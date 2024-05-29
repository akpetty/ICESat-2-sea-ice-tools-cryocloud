# coding: utf-8

"""
find_matching_IS2_S2_paris.py
Written by Marco Bagnardi (04/2022)

This script searches for semi-coincident ICESat-2 and Sentinel-2 data.
It uses the Sentinel API.

For each ICESat-2 granule in a given directory (e.g., data archive on 
cooler), the sofware find Sentinel-2 images with ICESat-2 data falling
within the image's footprint.

Output:
- PNG image showing the ICESat-2 granule name, the Sentinel-2 image
(or images) outline colorcoded by product name, the start and end time of
the data overlap, and the ICESat-2 footprint in black.
- TXT file with data summary and links for download of quicklook images or
full data products. Note that some of the Sentinel-2 may not be available
for direct download.
- Directory with all quicklook images of the Sentinel-2 tiles


UPDATE HISTORY:
Written 04/2022
Updated in 2023 by A Petty

Use the 

"""


from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
import datetime as dt
import h5py
import pandas as pd
import geopandas as gpd
import glob
from shapely.geometry import Point
import matplotlib.pyplot as plt
import os
from astropy.time import Time
import numpy as np
from PIL import Image
import io
import requests


# Function to download and save quicklook image to directory with IS-2 granule name
def save_quicklook(id_, username, password, dir_name, imagename):
    url = "https://apihub.copernicus.eu/apihub/odata/v1/Products('{}')/Products('Quicklook')/$value".format(id_)
    bytes_img = requests.session().get(url, auth=(username, password)).content
    ql = Image.open(io.BytesIO(bytes_img))

    ql = ql.save(dir_name + '/' + imagename + '.jpg')
    return


# Credentials for ESA catalog search
username = 'akpetty'
password = 'Icebridge01!'

#sentinelsat.exceptions.UnauthorizedError: Invalid user name or password. Note that account creation and password changes may take up to a week to propagate to the 'https://apihub.copernicus.eu/apihub/' API URL you are using. Consider switching to 'https://scihub.copernicus.eu/dhus/' instead in the mean time.
# Try https://apihub.copernicus.eu/apihub (Marco) or https://scihub.copernicus.eu/dhus/
# Initialize access to API
api = SentinelAPI(username, password, 'https://scihub.copernicus.eu/dhus/')

deltatime=20.0
maxcloud=10
lower_lat=60

# S2MSI1C: Level-1C data, is both radiometrically and geometrically corrected. This means that the images have been orthorectified, ensuring that they are properly aligned with real-world coordinates and free from distortions due to Earth's curvature and satellite perspective.
# S2MSI2A: Level-2A data, is an atmospherically corrected version of L1C data. It provides bottom-of-atmosphere (BOA) reflectance values, which have been adjusted for the effects of atmospheric gases and aerosols.
# Seems like most folk use S2MSI1C, maybe don't trust the atmospheric corrections??

sentinel2_product='S2MSI1C'

# ICESat-2 data location
ATL_input_dir = '/Users/aapetty/Data/ICESat2/r005/ATL07/'
ATL_filelist = glob.glob(ATL_input_dir + 'ATL07-01_20211*.h5')

# Loop over ICESat-2 granules
for ATL_filename in ATL_filelist:
    print(ATL_filename)

    ATL = h5py.File(ATL_filename, 'r')

    # Check spacecraft orientation and assign string beam IDs
    orientation = ATL['/orbit_info/sc_orient'][0]

    # Only use central strong beam locations
    if orientation == 0:
        beamID = 'gt2l'
    elif orientation == 1:
        beamID = 'gt2r'
    else:
        print('Spacecraft orientation not found.')

    # Extract data info from granule
    ATL_start_time = ATL['/ancillary_data/data_start_utc'][0]
    ATL_start = pd.to_datetime(ATL_start_time[:-8].decode('utf-8'), format='%Y-%m-%dT%H:%M:%S')
    ATL_end_time = ATL['/ancillary_data/data_end_utc'][0]
    ATL_end = pd.to_datetime(ATL_end_time[:-8].decode('utf-8'), format='%Y-%m-%dT%H:%M:%S')

    GPS_epoch = ATL['ancillary_data/atlas_sdp_gps_epoch'][:]

    # Build dataframe with location of data
    ATL_dF = pd.DataFrame({'Longitude': ATL[beamID + '/sea_ice_segments/longitude'][::200],
                           'Latitude': ATL[beamID + '/sea_ice_segments/latitude'][::200],
                           })
    ATL_dF['coords'] = list(zip(ATL_dF['Longitude'], ATL_dF['Latitude']))
    ATL_dF['coords'] = ATL_dF['coords'].apply(Point)

    GPS_time = ATL[beamID + '/sea_ice_segments/delta_time'][::200] + GPS_epoch
    
    # Use astropy to convert from gps time to datetime
    ATL_tgps = Time(GPS_time, format='gps')
    ATL_utc = ATL_tgps.utc.iso
    ATL_dF['UTC'] = ATL_utc

    ATL_gfd = gpd.GeoDataFrame(ATL_dF, geometry='coords')
    ATL_gfd = ATL_gfd.set_crs(4326, allow_override=True)

    
    #try:
    # Search for Sentinel-2 coincident data
    
    S2_query = api.query(platformname='Sentinel-2', producttype=sentinel2_product,
                         date=(ATL_start - dt.timedelta(minutes=deltatime), ATL_end + dt.timedelta(minutes=deltatime)),
                         cloudcoverpercentage=(0, maxcloud))

    

    S2_gdf = api.to_geodataframe(S2_query)

    S2_gdf_subset = S2_gdf[(S2_gdf.bounds.miny > lower_lat)]

    points_in_poly = gpd.tools.sjoin(ATL_gfd, S2_gdf_subset)
    #print(points_in_poly)
    
    try:
        # Empty geodataframes threw an exception here
        print('Number of overlapping tiles:', len(points_in_poly['title'].unique()))
    except:
        continue

    if len(points_in_poly['title'].unique()) > 0:

        # Filter products based on the tile ID
        #filtered_products = {k: v for k, v in products.items() if v['tileid'] == tile_id}

        print('all data:', S2_query)
        print('subset:', S2_gdf_subset)

        #print(S2_gdf_subset.title)

        # download all results from the search
        #api.download_all(S2_gdf_subset)

        

        #cwd = os.getcwd()
        #print('cwd:', cwd)

        atl_name = os.path.basename(ATL_filename)[:-3]
        filename = 'S2pairs_'+atl_name
        
        save_path = '/Users/aapetty/Data/S2/IS2_match/'+sentinel2_product+'/maxcloud'+str(maxcloud)+'_deltatime'+str(int(deltatime))+'/'+filename
        
        print('save_path:', save_path)
        
        if not os.path.exists(cwd+'/'+sentinel2_product):
             os.mkdir(cwd+'/'+sentinel2_product)

        if not os.path.exists(cwd+'/'+sentinel2_product+'/maxcloud'+str(maxcloud)+'_deltatime'+str(int(deltatime))):
             os.mkdir(cwd+'/'+sentinel2_product+'/maxcloud'+str(maxcloud)+'_deltatime'+str(int(deltatime)))

        if not os.path.exists(save_path):
             os.mkdir(save_path)

        
        f = open(save_path + "/" + filename+".txt", 'a')

        fig, ax = plt.subplots(1, 1)

        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(points_in_poly['title'].unique()))))
        y_shift = 0.8

        for granule in points_in_poly['title'].unique():
            print(granule)
            c = next(color)
            y_shift = y_shift -0.05

            S2_gdf_subset[S2_gdf_subset['title'] == granule].boundary.plot(ax=ax, color=c)
            plt.text(1.2, y_shift, str(granule), transform=ax.transAxes, color=c)

        points_in_poly.plot(ax=ax, markersize=3, color='k', marker='o')
        plt.title(os.path.basename(ATL_filename))
        plt.xlabel('Longitude (deg.)')
        plt.ylabel('Latitude (deg.)')

        plt.text(1.2, y_shift -0.10, points_in_poly.iloc[0].UTC, transform=ax.transAxes, color='k')
        plt.text(1.2, y_shift -0.15, points_in_poly.iloc[-1].UTC, transform=ax.transAxes, color='k')

        plt.savefig(save_path + "/" + filename+".png", bbox_inches='tight', dpi=100)

        f.write(save_path + '\n')
        f.write('\n')
        f.write('Lat min: ' + str(min(points_in_poly.Latitude)) + '\n')
        f.write('Lat max: ' + str(max(points_in_poly.Latitude)) + '\n')
        f.write('Time start: ' + points_in_poly.iloc[0].UTC + '\n')
        f.write('Time end: ' + points_in_poly.iloc[-1].UTC + '\n')
        f.write('\n')
        f.write(str(points_in_poly['title'].unique()) + '\n')
        f.write('\n')
        f.write(str(points_in_poly['summary'].unique()) + '\n')
        f.write('\n')
        f.write(str(points_in_poly['link'].unique()) + '\n')
        f.write('\n')
        f.write(str(points_in_poly['link_icon'].unique()) + '\n')
        f.write('\n')

        f.close()

        for title, uuid in zip(points_in_poly['title'].unique(), points_in_poly['uuid'].unique()):
            #print('saving quiklook tile / uuid:', title, uuid)
            #save_quicklook(uuid, username, password, save_path, title)

            print('Try downloading quicklook scenes using the API')
            try:
                api.download_quicklook(uuid, directory_path=save_path)
            except:
                print('Error downloading quicklook')
                continue

            print('Try downloading full scenes using the API')
            product_info = api.get_product_odata(uuid)
 
            is_online = api.is_online(uuid)

            if is_online:
                print(f'Product {uuid} is online. Starting download.')
                try:
                    # Added checksum=False as a weird exception related to the metadata
                    # https://github.com/sentinelsat/sentinelsat/issues/467
                    api.download(uuid, directory_path=save_path, checksum=False)
                except:
                    continue
            else:
                print(f'Product {uuid} is not online.')
                try:
                    api.trigger_offline_retrieval(uuid)
                except:
                    continue    

            
            #api.download(uuid, directory_path=save_path)

    #except:
    #    print('Some issue with the Sentinel-2 catalog has occurred, skipping current granule.')
