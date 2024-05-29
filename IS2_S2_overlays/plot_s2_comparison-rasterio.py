""" plot_s2_comparison.py
    
    Code to read in a Senintel-2 image, grab the ATL07/10 data and overaly/plot/analyze.
    Initial code written by Marco Bagnardi. Adapted by Alek Petty
  
    Output: 
      figures/release/profiles/S-2_overlay... - the main S-2/ICESat-2 overlay
      figures/release/profiles/s2-extra/croppedimage...: image rotated and cropped to be used by the profile overlay figure above
      figures/release/profiles/s2-extra/diagnostic...: looking at quality flag, adjusted freeboard and other metrics)
      figures/release/profiles/s2-extra/diagnostic_boxplot...: looking at number of height segmetns not included)
      figures/release/profiles/s2-extra/overlay...: the complete S-2 image with IS-2 overlay
      figures/release/profiles/s2-extra/AOI...: Area of Interest S-2 images

    example: 
      run as: python plot_s2_comparison.py -b 'gt1l' -e 'P5'

    Python dependencies:
        See below for the relevant module imports. More information on installation is given in the README file.

    Update history:
        08/30/2021: Version 1.
"""


### Import necessary Python modules
#import cartopy.crs as ccrs
#import cv2
try:
  from osgeo import gdal
  from osgeo import ogr
except:
  import gdal
  import ogr
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
from pyproj import Proj, transform
#import osr
from scipy import ndimage, misc
import sys
from glob import glob
import os
from matplotlib.ticker import MaxNLocator
import numpy.ma as ma
import argparse
import rasterio

# Our own module of bespoke functions
import utils as ut

def pick_profile(example):
  """ Pick S-2/ICESat-2 profiles to analyze"""

  if example=='P1':
  
    date='20190317_SH'

    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/RT_T20CPE_20190317T122209_B04.tif'
    # ATLAS DATA
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-02/data/ATL07-02_20190317*_12070201_'+relStr[-3:]+'_'+version+'.h5')[0]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-02/data/ATL10-02_20190317*12070201_'+relStr[-3:]+'_'+version+'.h5')[0]
    

    # Latitudinal limits
    Slim = -73
    Nlim = -72.57

    print('Arctic profile:', example, date, str(Slim)+'-'+str(Nlim))

    bufferRatio=60 # Bigger number smaler width (try 100 for zoom in)

  elif (example=='P2')|(example=='P3')|(example=='P4'):
      
      date='20190526_NH'

      version='*'
      
      # Sentinel-2 surface reflectance single-band geoTIFF file
      imageFile = s2_data_path+'/RT_T14XMQ_20190525T230121_B04.tif'
      # ATLAS DATA
      # -1 indicates we're using the latest version of the data (i.e. the higher number)
      ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190526*_08820301_'+relStr[-3:]+'_'+version+'.h5')[-1]
      ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190526*_08820301_'+relStr[-3:]+'_'+version+'.h5')[-1]

      if example=='P2':
        Slim = 80.57
        Nlim = 81

      if example=='P3':
        Slim = 80.15
        Nlim = 80.57

      if example=='P4':
        # focussing on the big lead in gt21 (need to also manually change the size of the top line when doing so!)
        Slim = 80.6
        Nlim = 80.66
      
      print('Arctic profile:', example, date, str(Slim)+'-'+str(Nlim))

      bufferRatio=60 # Bigger number = smaler width (try 100 for zoom in)


  elif (example=='P5'):
    #20190726_T09XWJ_T09XWH_T09XVH_T10XDQ_RGT441_DT22_AR

    date='20190726_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190726_T09XVH.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]

    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic summer profile:', example, date, str(Slim)+'-'+str(Nlim))


  elif (example=='P6'):
      
    date='20190726_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190726_T09XWJ.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    
    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic summer profile:', example, date, str(Slim)+'-'+str(Nlim))

  elif (example=='P7'):
    
    date='20190726_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190726_T09XWH.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    

    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic summer profile:', example, date, str(Slim)+'-'+str(Nlim))

  elif (example=='P8'):
    
    date='20190726_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190726_T10XDQ.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    
    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic summer profile:', example, date, str(Slim)+'-'+str(Nlim))

  elif (example=='P9'):

    date='20190622_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190622_T43XEK.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*1301*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*1301*'+relStr[-3:]+'_'+version+'.h5')[-1]
    
    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic late spring profile:', example, date, str(Slim)+'-'+str(Nlim))

  elif (example=='P10'):
      
    date='20190622_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190622_T43XEK.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*441*'+relStr[-3:]+'_'+version+'.h5')[-1]
    
    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic late spring profile:', example, date, str(Slim)+'-'+str(Nlim))

  elif (example=='P11'):
    #20190621_T11XNK_RGT1294_DT22_AR


    date='20190621_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/20190621_T11XNK.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/ATL07-01_20190726*1294*'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/ATL10-01_20190726*1294*'+relStr[-3:]+'_'+version+'.h5')[-1]
    
    # Just use very big bounds for now, should capture the entire image!
    Slim = 60
    Nlim = 88
    bufferRatio=80 # Bigger number smaler width (try 100 for zoom in)
    
    print('Arctic late spring profile:', example, date, str(Slim)+'-'+str(Nlim))


  #try:
  #  geoTIFF = rasterio.open(imageFile)
  #except IOError:
  #  sys.exit('geoTIFF file is not a valid file')

  ### Check if ICESat-2 ATLAS file is valid
  try:
    ATL07file = h5py.File(ATL07_dataFile, 'r')
  except IOError:
    sys.exit('not a valid file')

  ### Check if ICESat-2 ATLAS file is valid
  try:
    ATL10file = h5py.File(ATL10_dataFile, 'r')
  except IOError:
    sys.exit('not a valid file')


  return imageFile, ATL07file, ATL10file, Slim, Nlim, bufferRatio, date

def parse():
    """ Parse command line input arguments. """
    parser = argparse.ArgumentParser(description='S2 plot arguments')
    parser.add_argument('-e', action='store', default='', dest='example',
                        help='P1, P2, P3 in the paper, P4 is the zoomed in lead',
                        type=str)
    parser.add_argument('-b', action='store', default='', dest='beam_name',
                        help='beam name (e.g. gt1r)',
                        type=str)
    inps = parser.parse_args()
    return inps

def get_geoTIFF_info_rasterio(imageFile, plot=False):
  
  try:
    geoTIFF = rasterio.open(imageFile)
  except IOError:
    sys.exit('geoTIFF file is not a valid file')

  # Convert geoTIFF image to data array
  surfReflectance = geoTIFF.read()
  if surfReflectance.ndim>2:
    print('multi-band image so just take the 1st (red) band for now - WORK ON CONVERTING TO RGB INSTEAD!')
    surfReflectance=surfReflectance[0]
    band='Red'
  else:
    band='RGB'

  # Get extent of image (results slightly different to GDAL)
  geoTIFF_extent = geoTIFF.bounds

  # Extract spatial metadata
  geoTIFF_proj = geoTIFF.crs
  geoTIFF_geo_trans  = geoTIFF.transform

  #out_proj = 'epsg:' + str(geoTIFF_proj) # String with geoTIFF projection epsg code

  
  # Print basic GEOTIFF information
  print('Image size: ', geoTIFF.width, geoTIFF.height)
  print('Image projection: ', geoTIFF_proj)
  print('Xmin: ', geoTIFF_extent[0])
  print('Xmax: ', geoTIFF_extent[1])
  print('Ymin: ', geoTIFF_extent[2])
  print('Ymax: ', geoTIFF_extent[3])
  print('geoTIFF_extent', geoTIFF_extent)
  
  print('geoTIFF_geo_trans', geoTIFF_geo_trans)


  if plot==True:
    subplot_kw = dict(projection=geoTIFF_projection)
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot image as grayscale
    im = ax.imshow(surfReflectance, extent=geoTIFF_extent, origin='upper', cmap='gray')


  # Close dataset
  geoTIFF = None

  return surfReflectance, geoTIFF_extent, geoTIFF_proj, geoTIFF_geo_trans, band


def get_value_at_point(array_from, pos, geoTIFF_geo_trans):
    """ Extract raster value at given position from coordinates """
    try:
      samp_x = int((pos[0] - geoTIFF_geo_trans[2]) / geoTIFF_geo_trans[0])
      samp_y = int((pos[1] - geoTIFF_geo_trans[5]) / geoTIFF_geo_trans[4])

      return array_from[samp_y-1, samp_x-1]
    except:
      print('outside image extent', pos)

def createShapefile(filename, LL_x, LL_y, UL_x, UL_y, UR_x, UR_y, LR_x, LR_y):
    

    driver = ogr.GetDriverByName('ESRI Shapefile')

    datasource = driver.CreateDataSource(filename)
    layer = datasource.CreateLayer('layerName',geom_type=ogr.wkbPolygon)

    outline = ogr.Geometry(type=ogr.wkbLinearRing)
    outline.AddPoint(LL_x, LL_y)
    outline.AddPoint(UL_x, UL_y)
    outline.AddPoint(UR_x, UR_y)
    outline.AddPoint(LR_x, LR_y)
    outline.AddPoint(LL_x, LL_y)
    polygon = ogr.Geometry(type=ogr.wkbPolygon)
    polygon.AddGeometry(outline)

    #create feature object with polygon geometry type from layer object:
    feature = ogr.Feature( layer.GetLayerDefn() )
    feature.SetGeometry(polygon)
    layer.CreateFeature(feature)
   
    return

def rotate_IS2_S2(sub_dF07, bufferRatio):
  # Trigonometry calculations
  a = sub_dF07.UTM_Y.iloc[0] - sub_dF07.UTM_Y.iloc[-1]
  b = sub_dF07.UTM_X.iloc[0] - sub_dF07.UTM_X.iloc[-1]
  c = np.sqrt((a**2) + (b**2))

  cosA = (b**2 + c**2 - a**2) / (2*b*c)
  inclin = np.arccos(cosA)
  inclin_deg = np.degrees(inclin)

  # Adjust size of buffer zone from ICESat-2 ground track depending on length of subset
  buffer = c/bufferRatio

  deltaX = buffer * np.sin(inclin)
  deltaY = buffer * np.cos(inclin)

  # Compute coordinates of corners of AOI

  # Descending satellite pass
  if sub_dF07.UTM_Y.iloc[0] > sub_dF07.UTM_Y.iloc[-1]:
    rev = True # set flag to determine how to apply rotation
    
    UL_x = sub_dF07.UTM_X.iloc[0] - deltaX
    UL_y = sub_dF07.UTM_Y.iloc[0] + deltaY

    UR_x = sub_dF07.UTM_X.iloc[0] + deltaX
    UR_y = sub_dF07.UTM_Y.iloc[0] - deltaY

    LL_x = sub_dF07.UTM_X.iloc[-1] - deltaX
    LL_y = sub_dF07.UTM_Y.iloc[-1] + deltaY

    LR_x = sub_dF07.UTM_X.iloc[-1] + deltaX
    LR_y = sub_dF07.UTM_Y.iloc[-1] - deltaY

  # Ascending satellite pass    
  if sub_dF07.UTM_Y.iloc[0] < sub_dF07.UTM_Y.iloc[-1]:
    rev = False # set flag to determine how to apply rotation
    
    UL_x = sub_dF07.UTM_X.iloc[0] - deltaX
    UL_y = sub_dF07.UTM_Y.iloc[0] - deltaY

    UR_x = sub_dF07.UTM_X.iloc[0] + deltaX
    UR_y = sub_dF07.UTM_Y.iloc[0] + deltaY

    LL_x = sub_dF07.UTM_X.iloc[-1] - deltaX
    LL_y = sub_dF07.UTM_Y.iloc[-1] - deltaY

    LR_x = sub_dF07.UTM_X.iloc[-1] + deltaX
    LR_y = sub_dF07.UTM_Y.iloc[-1] + deltaY

  # Print AOI information
  print('ICESat-2 inclination angle in imagery coordinate system: ', inclin_deg)
  print('AOI coordinates:')
  print('UL: ', UL_x, UL_y)
  print('UR: ', UR_x, UR_y)
  print('LR: ', LR_x, LR_y)
  print('LL: ', LL_x, LL_y)

  return rev, inclin_deg, UL_x, UL_y, UR_x, UR_y, LR_x, LR_y, LL_x, LL_y, buffer

def plot_height_freeboard(sub_dF07, sub_dF10):

  # play around with adjusting the bounds on the freeboard data (removing leads etc)
  sub_dF10_adj = sub_dF10[(sub_dF10['ssh_flag'] < 0.5) & (sub_dF10['seg_type_flag']==1 ) & (sub_dF10['freeboard']>0.04)]
  
  mean_fb_adj=np.round(np.mean(sub_dF10_adj['freeboard']), decimals=3)

  shadingSize=30
  cmap_ssh = colors.ListedColormap(['red', 'yellow'])

  # Generate plot with white background
  #%matplotlib qt
  #%matplotlib inline
  fig= plt.figure(figsize=(10,8))
  #fig.patch.set_facecolor('xkcd:white')

  # ATL07 heights ssh label
  ax1 = plt.subplot2grid((5, 1), (0, 0), rowspan=1)

  for index, row in sub_dF07.iterrows():
      #x0 = row['along_dist'] - row['seg_length']/2000
      x0 = row['along_dist'] - shadingSize/2000
      x1 = row['along_dist'] + shadingSize/2000
      
      if row['ssh_flag'] > 0.1:
          plt.axvspan(x0,x1, facecolor='k', alpha=0.2)
          
  plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=sub_dF07['height'], cmap='hot', s=2, zorder=2)
  plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
  ax1.xaxis.set_ticklabels([])
  plt.ylim(-0.5, 3)
  ax1.annotate('Grey shading is ssh_flag',xy=(0.02, 0.87), xycoords='axes fraction')
  plt.ylabel('ATL07 Height')
  # ATL07 heights (in atl10)
  ax2 = plt.subplot2grid((5, 1), (1, 0), rowspan=1)

  for index, row in sub_dF07.iterrows():
      x0 = row['along_dist'] - shadingSize/2000
      x1 = row['along_dist'] + shadingSize/2000
      
      if row['ssh_flag'] > 0.1:
          plt.axvspan(x0,x1, facecolor='k', alpha=0.2)
          
  plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=atlmaskint, vmin=0, vmax=1, cmap='viridis_r', s=2, zorder=2)
  plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
  plt.ylim(-0.5, 3)
  ax2.xaxis.set_ticklabels([])
  plt.ylabel('ATL07 Height')
  ax2.annotate('Purple markers not in ATL10',xy=(0.02, 0.87), xycoords='axes fraction')

  ax3 = plt.subplot2grid((5, 1), (2, 0), rowspan=1)
  #for x in range(np.size(x0)):
      #print(x)
  #for m in range(np.size(atlmask)):
  #    if atlmask[m]==False:
  #        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)

  for index, row in sub_dF07.iterrows():
      x0 = row['along_dist'] - shadingSize/2000
      x1 = row['along_dist'] + shadingSize/2000
      
      if row['ssh_flag'] > 0.1:
          plt.axvspan(x0,x1, facecolor='k', alpha=0.2)


  cmap = plt.cm.get_cmap('gnuplot', 7)    # 11 discrete colors

  imq=plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=sub_dF07['quality_flag'], cmap=cmap, vmin=-1.5, vmax=5.5, s=2, zorder=2)
  plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
  plt.ylim(-0.5, 3)
  ax3.xaxis.set_ticklabels([])
  plt.ylabel('ATL07 Height')
  percent=np.round(np.size(np.where(sub_dF07['quality_flag']>4))/np.size(sub_dF07['quality_flag'])*100, decimals=2)

  ax3.annotate('Quality flag ('+str(percent)+'% > 4)',xy=(0.02, 0.87), xycoords='axes fraction')

  cax = fig.add_axes([0.925, 0.41, 0.022, 0.2])
  cbar = plt.colorbar(imq,cax=cax, orientation='vertical', extend='both', use_gridspec=True)
  #cbar.set_label('Quality flag', labelpad=28, rotation=0)
  xticks1 = np.linspace(-1, 5, 7)
  cbar.set_ticks(xticks1)
  cbar.set_ticklabels(['inv', '0', '1', '2', '3', '4', '5'])
  #plt.clim(-1.5,5.5)

  ax4 = plt.subplot2grid((5, 1), (3, 0), rowspan=1)
  #mask = np.isin(sub_DF['height_segment_id'].values, sub_DF10['height_segment_id'].values,  invert=True)
  #for m in range(np.size(atlmask)):
  #    if atlmask[m]==False:
  #        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)
  # Add lead_flag=1 shading
  for index, row in sub_dF10.iterrows():
      # CHANGE TO LEAD LENGTH
      x0 = row['along_dist'] - shadingSize/2000
      x1 = row['along_dist'] + shadingSize/2000
      
      if row['ssh_flag'] > 0.1:
          plt.axvspan(x0,x1, facecolor='k', alpha=0.2)

  ax4.annotate('Grey shading is leads',xy=(0.02, 0.87), xycoords='axes fraction')
  plt.scatter(sub_dF10['along_dist'], sub_dF10['freeboard'], c=sub_dF10['freeboard'], cmap='hot', s=2, zorder=2)
  #ax3.xaxis.set_ticklabels([])
  plt.ylabel('Freeboard')
  plt.xlabel('Along track distance (km)')
  plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
  plt.ylim(-0.5, 3)
  ax4.annotate('Mean = '+str(mean_fb)+' m',xy=(0.8, 0.84), xycoords='axes fraction')

  ax5 = plt.subplot2grid((5, 1), (4, 0), rowspan=1)
  #mask = np.isin(sub_DF['height_segment_id'].values, sub_DF10['height_segment_id'].values,  invert=True)
  #for m in range(np.size(atlmask)):
  #    if atlmask[m]==False:
  #        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)
  # Add lead_flag=1 shading
  for index, row in sub_dF10.iterrows():
      # CHANGE TO LEAD LENGTH
      x0 = row['along_dist'] - shadingSize/2000
      x1 = row['along_dist'] + shadingSize/2000
      
      if row['ssh_flag'] > 0.1:
          plt.axvspan(x0,x1, facecolor='k', alpha=0.2)

  ax5.annotate('Adjusted freeboard',xy=(0.02, 0.87), xycoords='axes fraction')
  plt.scatter(sub_dF10_adj['along_dist'], sub_dF10_adj['freeboard'], c=sub_dF10_adj['freeboard'], cmap='hot', s=2, zorder=2)
  #ax3.xaxis.set_ticklabels([])
  plt.ylabel('Freeboard')
  plt.xlabel('Along track distance (km)')
  plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
  plt.ylim(-0.5, 3)
  ax5.annotate('Mean = '+str(mean_fb_adj)+' m',xy=(0.8, 0.84), xycoords='axes fraction')


  plt.subplots_adjust(left=0.06, right=0.92, top=0.96, bottom=0.08, wspace=0, hspace=0)
  #plt.tight_layout()
  ### Save plot to file
  plt.savefig(figPath+'s2-extra/diagnostic_'+date+beam_name+relStr+str(Slim)+str(Nlim)+example+'.png', dpi=300)

def plot_diagnostic_boxplot(sub_dF07):
  f, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (11, 4), sharex=False, sharey=True)
  plt.sca(ax1)
  fig.patch.set_facecolor('xkcd:white')
  plt.scatter(sub_dF07['along_dist'].values[atlmask], sub_dF07['height'].values[atlmask], color='r', s=2, alpha=0.3)
  #plt.scatter(sub_DF07['along_dist'].values[~atlmask], sub_DF07['elev'].values[~atlmask], color='b', s=2, alpha=0.3)
  plt.annotate('Segments not in ATL10', xy=(0.02, 0.92), xycoords='axes fraction')
  plt.xlabel('Along track distance (km)')
  plt.ylabel('Heights (m)')
  plt.sca(ax2)
  plt.boxplot([sub_dF07['height'].values[~np.isnan(sub_dF07['height'].values)], sub_dF07['height'].values[~atlmask][~np.isnan(sub_dF07['height'].values[~atlmask])], sub_dF07['height'].values[atlmask][~np.isnan(sub_dF07['height'].values[atlmask])]])
  ax2.set_xticklabels(['ATL07 Heights', 'In ATL10', 'Not in ATL10'])
  plt.savefig(figPath+'s2-extra/diagnostic_boxplot_'+date+beam_name+relStr+str(Slim)+str(Nlim)+example+'.png', dpi=300)

### 1. Data ingestion and preparation
################### Input data file names ########################################

relStr='rel006'
data_path = '/home/jovyan/Data/IS2/'

s2_data_path = '/home/jovyan/Data/S2/'

figPath='/home/jovyan/GitHub-output/ICESat-2-sea-ice-tools/plot_s2_comparison'

if not os.path.exists(figPath):
  os.makedirs(figPath)
if not os.path.exists(figPath+'s2-extra/'):
  os.makedirs(figPath+'s2-extra/')
#if not os.path.exists(savePath):
#  os.makedirs(savePath)

################### Read command line input parameters ########################################

inps = parse()
# Check for missing command line input parameters
if inps.example is None:
    print('Missing example profile index. Use option -e')

    sys.exit(1)
if inps.beam_name is None:
    print('Missing beam. Use option -b')
    sys.exit(1)

example = inps.example
beam_name = inps.beam_name

print(example, beam_name)

band=0 #0 = red (in some cases set to np.nan if just single band)

################### Read S-2/ICESat-2 data for given example profile ########################################
print('Read in S-2/ICESat-2 data...')
imageFile, ATL07_dataFile, ATL10_dataFile, Slim, Nlim, bufferRatio, date = pick_profile(example)

### Read geoTIFF file and print spatial extent information
surfReflectance, geoTIFF_extent, out_proj, geoTIFF_geo_trans, band = get_geoTIFF_info_rasterio(imageFile)

### Read ATL07 file into pandas dataframe
ATL07dF, beamStr, beamStrength = ut.get_atl07_data_beam_extra(ATL07_dataFile, beamStr=beam_name, km_unit=True)
ATL10dF, _ = ut.get_atl10_data_beam_extra(ATL10_dataFile, beamStr=beam_name, km_unit=True)

### Transform WGS84 Latitude and Longitude to UTM X and Y
ATL_inProj = Proj(init='epsg:4326') # ATLAS data projection
ATL_outProj = Proj(init=out_proj) # GeotTIFF projection

seg_x_utm, seg_y_utm = transform(ATL_inProj, ATL_outProj, ATL07dF['lons'].values, ATL07dF['lats'].values) # Transform coordinates
seg_x_utm10, seg_y_utm10 = transform(ATL_inProj, ATL_outProj, ATL10dF['lons'].values, ATL10dF['lats'].values) # Transform coordinates

# Add new coorinates to dataframe
ATL07dF['UTM_X'] = seg_x_utm
ATL07dF['UTM_Y'] = seg_y_utm

# Add new coorinates to dataframe
ATL10dF['UTM_X'] = seg_x_utm10
ATL10dF['UTM_Y'] = seg_y_utm10

# Replace no-data value 3.4028235e+38 with Numpy NaN in data frame
ATL07dF = ATL07dF.replace(np.max(ATL07dF.height), np.nan)
# Replace no-data value 3.4028235e+38 with Numpy NaN in data frame
#ATL10dF = ATL10dF.replace(np.max(ATL10dF.freeboard), np.nan)

# Crop ICESat-2 data frame to imagery extent
ATL07dF_crop = ATL07dF[(ATL07dF['UTM_X'] > geoTIFF_extent[0]) & (ATL07dF['UTM_X'] < geoTIFF_extent[1]) & 
                   (ATL07dF['UTM_Y'] > geoTIFF_extent[2]) & (ATL07dF['UTM_Y'] < geoTIFF_extent[3])]
# Crop ICESat-2 data frame to imagery extent
ATL10dF_crop = ATL10dF[(ATL10dF['UTM_X'] > geoTIFF_extent[0]) & (ATL10dF['UTM_X'] < geoTIFF_extent[1]) & 
                   (ATL10dF['UTM_Y'] > geoTIFF_extent[2]) & (ATL10dF['UTM_Y'] < geoTIFF_extent[3])]

# Reset data frame index value
ATL07dF_crop = ATL07dF_crop.reset_index(drop=True)
# Reset data frame index value
ATL10dF_crop = ATL10dF_crop.reset_index(drop=True)

# Print dataframe size 
print('Full ATL07 dataset n. of columns: ', ATL07dF.shape[1])
print('Full ATL07 dataset n. of segments: ', ATL07dF.shape[0])
print('Cropped ATL07 dataset n. of columns: ', ATL07dF_crop.shape[1])
print('Cropped ATL07 dataset n. of segments: ', ATL07dF_crop.shape[0])

# Print ATL10 dataframe size 
print('Full ATL10 dataset n. of columns: ', ATL10dF.shape[1])
print('Full ATL10 dataset n. of segments: ', ATL10dF.shape[0])
print('Cropped ATL10 dataset n. of columns: ', ATL10dF_crop.shape[1])
print('Cropped ATL10 dataset n. of segments: ', ATL10dF_crop.shape[0])

### Extract Sentinel-2 surface reflectance values at ICESat-2 segment locations
print('Extract S-2 data along ICESat-2 beam...')

# Preallocate empty array
dF_length = ATL07dF_crop.shape[0]
band_value = np.empty([dF_length, 1])

# Extract imagery pixel values at locations of ICESat-2 data
for i in range(dF_length):
    band_value[i] = get_value_at_point(surfReflectance, (ATL07dF_crop['UTM_X'][i], ATL07dF_crop['UTM_Y'][i]), geoTIFF_geo_trans)

# Add new column to data frame
ATL07dF_crop['BandValue'] = band_value


### Crop data to AOI
sub_dF07 = ATL07dF_crop[(ATL07dF_crop['lats'] < Nlim) & (ATL07dF_crop['lats'] > Slim) ]
# Reset data frame index value
sub_dF07 = sub_dF07.reset_index(drop=True)

# Print dataframe size 
print('AOI n. of columns: ', sub_dF07.shape[1], " * added surface reflectance value.")
print('AOI n. of points: ', sub_dF07.shape[0])


### Subset datasets to limits
sub_dF10 = ATL10dF_crop[(ATL10dF_crop['lats'] < Nlim) & (ATL10dF_crop['lats'] > Slim) ]

# Reset data frame index value
sub_dF10 = sub_dF10.reset_index(drop=True)

# Print dataframe size 
print('AOI n. of columns: ', sub_dF10.shape[1], " * added surface refelctance value.")
print('AOI n. of points: ', sub_dF10.shape[0])

print('start lat/lon of subset', sub_dF07['lats'].iloc[0], sub_dF07['lons'].iloc[0])
print('end lat/lon of subset', sub_dF07['lats'].iloc[-1], sub_dF07['lons'].iloc[-1])

print('start lat/lon of ATL10 subset', sub_dF10['lats'].iloc[0], sub_dF10['lons'].iloc[0])
print('end lat/lon of ATL10 subset', sub_dF10['lats'].iloc[-1], sub_dF10['lons'].iloc[-1])

print('start lat/lon of profile across entire image', ATL07dF_crop['lats'].iloc[0], ATL07dF_crop['lons'].iloc[0])
print('end lat/lon of profile across entire image', ATL07dF_crop['lats'].iloc[-1], ATL07dF_crop['lons'].iloc[-1])

print('length of subset', str(int(0.001*(sub_dF07['along_dist'].iloc[-1] - sub_dF07['along_dist'].iloc[0])))+' km')
print('length of ATL10 subset', str(int(0.001*(sub_dF10['along_dist'].iloc[-1] - sub_dF10['along_dist'].iloc[0])))+' km')
print('length of profile across entire image', str(int(0.001*(ATL07dF_crop['along_dist'].iloc[-1] - ATL07dF_crop['along_dist'].iloc[0])))+' km')


### Determine ICESat-2 ground track inclination in Sentinel-2 image

rev, inclin_deg, UL_x, UL_y, UR_x, UR_y, LR_x, LR_y, LL_x, LL_y, buffer  = rotate_IS2_S2(sub_dF07, bufferRatio)

### Make new plot of Sentinel-2 image with ICESat-2 beam overlaid (AOI only)

# Select ICESat-2 parameter to plot
IS2_param = 'ssh_flag'

#subplot_kw = dict(projection=geoTIFF_projection)
fig, ax = plt.subplots(figsize=(10, 10))

# Plot image as grayscale
im = ax.imshow(surfReflectance, extent=geoTIFF_extent, origin='upper', cmap='gray')

# Overlay ICESat-2 data
plt.scatter(sub_dF07['UTM_X'], sub_dF07['UTM_Y'], s=1, c=sub_dF07[IS2_param], cmap='viridis')
plt.plot([UL_x, LL_x, LR_x, UR_x, UL_x], [UL_y, LL_y, LR_y, UR_y, UL_y], c='r')
plt.colorbar()
plt.savefig(figPath+'s2-extra/overlay'+date+beam_name+relStr+str(Slim)+str(Nlim)+example+'.png', dpi=300)

### Create shapefile with extent of AOI to be used for image cropping
# Call function to create AOI polygon shapefile
createShapefile(figPath+'s2-extra/AOIPolygon.shp', LL_x, LL_y, UL_x, UL_y, UR_x, UR_y, LR_x, LR_y)
shapefile=figPath+'AOIPolygon.shp'
outpath=figPath+'AOI.tif'

### Use GDAL to crop Sentinel-2 image to AOI

ds = gdal.Warp(figPath+'s2-extra/AOI.tif',
               imageFile,
               format='GTiff',
               cutlineDSName= figPath+'s2-extra/AOIPolygon.shp',
               cropToCutline = True,
               dstAlpha = True)
ds = None

### Convert GEOTiff AOI file to JPG for subsequent use

# Options to convert GEOTiff image to JPEG
options_list = ['-ot Byte','-of JPEG','-b 1','-scale'] 
options_string = " ".join(options_list)

# Use gdal_translate in Python
gdal.Translate(figPath+'s2-extra/AOI.jpg',
               figPath+'s2-extra/AOI.tif',
               options=options_string)

# Read new JPEG image
#image = cv2.imread('AOI.tif')
image = plt.imread(figPath+'s2-extra/AOI.jpg')


# Rotate image by inclination angle according to satellite pass
# Descending
if rev is True:
    rot_im = ndimage.rotate(image, 180-inclin_deg, reshape=True, order=1)
# Ascending
if rev is False:
    rot_im = ndimage.rotate(image, inclin_deg, reshape=True, order=1)

### Trim edges to remove image padding
padNum=12
crop = np.delete(rot_im,np.where(~rot_im.any(axis=0))[0], axis=1)
crop = np.delete(crop,np.where(~crop.any(axis=1))[0], axis=0)
crop = np.delete(crop, range(padNum), 0)
crop = np.delete(crop, range(crop.shape[0]-padNum,crop.shape[0]), 0)
crop = np.delete(crop, range(padNum), 1)
crop = np.delete(crop, range(crop.shape[1]-padNum,crop.shape[1]), 1)

# Need to add a bigger loop to just load this if already generated the above!
#crop.dump(figPath+'s2-extra/cropped_image'+date+beam_name+relStr+str(Slim)+str(Nlim)+example)
np.savetxt(figPath+'s2-extra/cropped_image'+date+beam_name+relStr+str(Slim)+str(Nlim)+example+'.txt', crop)
print(crop.shape)


### Prepare data for plotting
# Filp image if necessary
if rev is False:
    crop = np.flipud(np.fliplr(crop))

start_index = sub_dF07['along_dist'].iloc[0]
# Convert distances from meters to kilometers with respect to first segment in AOI
sub_dF07['along_dist'] = sub_dF07['along_dist'] - start_index
#sub_dF07['along_dist'] = sub_dF07['along_dist'].multiply(0.001)

# Convert distances from meters to kilometers with respect to first segment in AOI
sub_dF10['along_dist'] = sub_dF10['along_dist'] - start_index
#sub_dF10['along_dist'] = sub_dF10['along_dist'].multiply(0.001)


# Generate and print some stats about this subset
atl07_segs=np.size(sub_dF07['ssh_flag'].values)
atl07_segs

ssh07_segs=np.size(np.where(sub_dF07['ssh_flag'] == 1))
ssh07_segs

specular_segs=np.size(np.where((sub_dF07['seg_type'] > 1.5)&(sub_dF07['seg_type'] < 5.5)))
specular_segs

darklead_segs=np.size(np.where(sub_dF07['seg_type'] > 5.5))
darklead_segs

cloud_segs=np.size(np.where(sub_dF07['seg_type'] < 0.5))
cloud_segs

atl10_segs=np.size(sub_dF10['ssh_flag'].values)
atl10_segs

ssh10_segs=np.size(np.where(sub_dF10['ssh_flag'] > 1.5))
ssh10_segs

mean_fb=np.round(np.mean(sub_dF10['freeboard']), decimals=3)
mean_fb

mafreeboards=ma.masked_where(sub_dF10['freeboard']<0.06, sub_dF10['freeboard'])

# Find the segments that are in ATL07 but not ATL10 
atlmask = np.isin(sub_dF07['height_id'].values, sub_dF10['height_segment_id'].values,  invert=True)
atlmaskint=atlmask.astype(int)
np.size(atlmask), np.size(np.nonzero(atlmaskint)[0])

#ssh ids in ATL10
ssh_ids=sub_dF07['height_id'].values[np.where(sub_dF07['ssh_flag'].values>0.5)]
ssh_inatl10 = np.isin(ssh_ids, sub_dF10['height_segment_id'].values)

#ssh ids
radio_ids=sub_dF07['height_id'].values[np.where(sub_dF07['seg_type'].values<1.5)]
radio_inatl10 = np.isin(radio_ids, sub_dF10['height_segment_id'].values,  invert=True)
radio_inatl10

percent=np.round(np.size(np.where(sub_dF07['quality_flag']==5))/np.size(sub_dF07['quality_flag'])*100, decimals=2)
percent

# Plot ATL07 diagnostic boxplot
plot_diagnostic_boxplot(sub_dF07)

# Plot the height/freeboard data
plot_height_freeboard(sub_dF07, sub_dF10)


print('plotting final profile...')
shadingSize=30


#=========== Subplot 1 ===========
fig= plt.figure(figsize=(9,11))
fig.patch.set_facecolor('xkcd:white')

ax1 = plt.subplot2grid((14, 1), (0, 0), rowspan=2)
plt.imshow(crop, extent=[sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1], -buffer/1000, +buffer/1000], aspect='auto', cmap='gray', vmin=0, vmax=255)
plt.scatter(sub_dF07['along_dist'][sub_dF07['ssh_flag']<0.5],np.zeros(len(sub_dF07['along_dist'][sub_dF07['ssh_flag']<0.5])), c='r', s=buffer/800)
plt.scatter(sub_dF07['along_dist'][sub_dF07['ssh_flag']>0.5],np.zeros(len(sub_dF07['along_dist'][sub_dF07['ssh_flag']>0.5])), c='y', s=buffer/400, zorder=2)

#plt.scatter(sub_dF07['along_dist'],np.zeros(len(sub_dF07['along_dist'])), c=sub_dF07['ssh_flag''], s=buffer/10000, cmap=cmap_ssh)

ax1.xaxis.set_ticklabels([])
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])
plt.autoscale(enable=True, axis='x', tight=True)
#plt.suptitle(date, color='k', fontsize=12)
ax1.annotate(date+' '+relStr+' '+beam_name+' ('+beamStrength+')  N(ATL07) segments = '+str(atl07_segs)+' N(ATL10) = '+str(atl10_segs),xy=(0.02, 1.05), xycoords='axes fraction')

#=========== Subplot 2 ===========
ax2 = plt.subplot2grid((14, 1), (2, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if (row['ssh_flag'] > 0.1):
        plt.axvspan(x0,x1, facecolor='k', alpha=0.1)

plt.scatter(sub_dF07['along_dist'],sub_dF07['BandValue'], c=sub_dF07['BandValue'], s=1, cmap='gray', vmin = 0, vmax = 255, zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax2.xaxis.set_ticklabels([])
ax2.set_yticks([0, 255])
plt.ylim(-1,266)
plt.ylabel('Red band', color='k', fontsize=12)
ax2.tick_params(axis='y', colors='k', labelsize=12)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])

#=========== Subplot 3 ===========
ax3 = plt.subplot2grid((14, 1), (4, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # Specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    # Dark lead
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['height'], c=sub_dF07['height'], s=1, vmin=0, vmax=3, cmap='hot', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax3.xaxis.set_ticklabels([])
ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
plt.ylabel('Height (m)', color='k', fontsize=12)
ax3.tick_params(axis='y', colors='k', labelsize=12)
plt.ylim(-0.2,4)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])

ax3.annotate('N='+str(specular_segs)+' (specular segs)',xy=(0.02, 0.87), color='b', xycoords='axes fraction')
ax3.annotate('N='+str(darklead_segs)+' (dark lead segs)',xy=(0.02, 0.75), color='r', xycoords='axes fraction')
ax3.annotate('N='+str(cloud_segs)+' (cloud segs)',xy=(0.02, 0.63), color='y', xycoords='axes fraction')
ax3.annotate('N='+str(ssh07_segs)+' (ATL07 candidate ssh segs)',xy=(0.02, 0.51), xycoords='axes fraction')

#=========== Subplot 4 ===========
ax4 = plt.subplot2grid((14, 1), (6, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    #x0 = row['along_dist'] - row['seg_length']/2000
    #x1 = row['along_dist'] + row['seg_length']/2000
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['photon_rate'], c=sub_dF07['photon_rate'], s=1, cmap='cool', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax4.xaxis.set_ticklabels([])

plt.ylabel('Photon Rate (ph/shot)', color='k', fontsize=12)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])

ax4.tick_params(axis='y', colors='k', labelsize=12)
ax4.yaxis.set_major_locator(MaxNLocator(5))

#=========== Subplot 5 ===========
ax5 = plt.subplot2grid((14, 1), (8, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['bground_rate'], c=sub_dF07['bground_rate'], s=1, cmap='summer', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax5.xaxis.set_ticklabels([])
ax5.yaxis.set_label_position("right")
ax5.yaxis.tick_right()
plt.ylabel('Background Rate (MHz)', color='k', fontsize=12)
plt.xlabel('Along track distance (km)', color='k', fontsize=16)
ax5.tick_params(axis='y', colors='k', labelsize=12)
ax5.tick_params(axis='x', colors='k', labelsize=12)
ax5.yaxis.set_major_locator(MaxNLocator(integer=True))

#=========== Subplot 6 ===========
ax6 = plt.subplot2grid((14, 1), (10, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # Specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    # Dark lead
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],100.*(sub_dF07['segs_used']/sub_dF07['segs_total']), c=100.*(sub_dF07['segs_used']/sub_dF07['segs_total']), s=1, vmin=0, vmax=100, cmap='spring', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax6.xaxis.set_ticklabels([])
#ax6.yaxis.set_label_position("left")
#ax6.yaxis.tick_left()
plt.ylabel('Pulses used (%)', color='k', fontsize=12)
ax6.tick_params(axis='y', colors='k', labelsize=12)
#plt.ylim(-10,4)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])


#=========== Subplot 7 ===========
ax7 = plt.subplot2grid((14, 1), (12, 0), rowspan=2)
           
plt.scatter(sub_dF07['along_dist'],sub_dF07['cloud_flag_asr'], c=sub_dF07['layer_flag'], s=2, vmin=-0.1, vmax=1.1, cmap='RdBu', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax7.xaxis.set_ticklabels([])
#ax6.yaxis.set_label_position("left")
#ax6.yaxis.tick_left()
plt.ylabel('ASR cloud flag (0=clear, 5=cloudy)', color='k', fontsize=10)
ax7.tick_params(axis='y', colors='k', labelsize=12)
ax7.yaxis.set_label_position("right")
ax7.yaxis.tick_right()
plt.ylim(-0.5,6.5)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])
#plt.xlabel('Along track distance (km)', color='k', fontsize=12)
#ax7.xaxis.set_ticklabels([])
ax7.annotate(r'colored by aggregate layer_flag (red=clear) or (blue=cloudy)',xy=(0.01, 0.88), color='k',xycoords='axes fraction')
ax7.xaxis.set_ticklabels([])

plt.xlabel('Along track distance (km)', color='k', fontsize=12)
#ax7.xaxis.set_ticklabels([])


plt.subplots_adjust(left=0.08, right=0.94, top=0.96, bottom=0.06, wspace=0, hspace=0)
### Save plot to file
plt.savefig(figPath+'/S-2_overlay_'+date+beam_name+relStr+str(Slim)+str(Nlim)+example+'.png', dpi=300)






#=========== Subplot 1 ===========
fig= plt.figure(figsize=(9,7))
fig.patch.set_facecolor('xkcd:white')

ax1 = plt.subplot2grid((10, 1), (0, 0), rowspan=2)
plt.imshow(crop, extent=[sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1], -buffer/1000, +buffer/1000], aspect='auto', cmap='gray', vmin=0, vmax=255)
plt.scatter(sub_dF07['along_dist'][sub_dF07['ssh_flag']<0.5],np.zeros(len(sub_dF07['along_dist'][sub_dF07['ssh_flag']<0.5])), c='r', s=buffer/800)
plt.scatter(sub_dF07['along_dist'][sub_dF07['ssh_flag']>0.5],np.zeros(len(sub_dF07['along_dist'][sub_dF07['ssh_flag']>0.5])), c='y', s=buffer/400, zorder=2)

#plt.scatter(sub_dF07['along_dist'],np.zeros(len(sub_dF07['along_dist'])), c=sub_dF07['ssh_flag''], s=buffer/10000, cmap=cmap_ssh)

ax1.xaxis.set_ticklabels([])
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])
plt.autoscale(enable=True, axis='x', tight=True)
#plt.suptitle(date, color='k', fontsize=12)
ax1.annotate(date+' '+relStr+' '+beam_name+' ('+beamStrength+')  N(ATL07) segments = '+str(atl07_segs)+' N(ATL10) = '+str(atl10_segs),xy=(0.02, 1.05), xycoords='axes fraction')

#=========== Subplot 2 ===========
ax2 = plt.subplot2grid((10, 1), (2, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if (row['ssh_flag'] > 0.1):
        plt.axvspan(x0,x1, facecolor='k', alpha=0.1)

plt.scatter(sub_dF07['along_dist'],sub_dF07['BandValue'], c=sub_dF07['BandValue'], s=1, cmap='gray', vmin = 0, vmax = 255, zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax2.xaxis.set_ticklabels([])
if band=='RGB':
  ax2.set_yticks([0, 1])
  plt.ylim(-0.1,1.1)
  plt.ylabel('RGB', color='k', fontsize=12)
else:
  ax2.set_yticks([0, 255])
  plt.ylim(-1,266)
  plt.ylabel('Red band', color='k', fontsize=12)

ax2.tick_params(axis='y', colors='k', labelsize=12)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])

#=========== Subplot 3 ===========
ax3 = plt.subplot2grid((10, 1), (4, 0), rowspan=2)


for index, row in sub_dF07.iterrows():
    #x0 = row['along_dist'] - row['seg_length']/2000
    #x1 = row['along_dist'] + row['seg_length']/2000
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['photon_rate'], c=sub_dF07['photon_rate'], s=1, cmap='hot', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax3.xaxis.set_ticklabels([])

plt.ylabel('Photon Rate (ph/shot)', color='k', fontsize=12)
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])
ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
ax3.tick_params(axis='y', colors='k', labelsize=12)
ax3.yaxis.set_major_locator(MaxNLocator(5))

#=========== Subplot 4 ===========
ax4 = plt.subplot2grid((10, 1), (6, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['bground_rate_norm'], c=sub_dF07['bground_rate_norm'], s=1, cmap='summer', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax4.xaxis.set_ticklabels([])

plt.ylabel('Norm Bground Rate (MHz)', color='k', fontsize=12)
#
ax4.tick_params(axis='y', colors='k', labelsize=12)
ax4.tick_params(axis='x', colors='k', labelsize=12)
ax4.yaxis.set_major_locator(MaxNLocator(integer=True))

#=========== Subplot 5 ===========
ax5 = plt.subplot2grid((10, 1), (8, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # Specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    # Dark lead
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)

width_max = np.nanpercentile(sub_dF07['height_width'].values, 99)
plt.scatter(sub_dF07['along_dist'],sub_dF07['height_width'], c=sub_dF07['height_width'], s=1, vmin=0, vmax=width_max, cmap='winter', zorder=2)
#plt.autoscale(enable=True, axis='x', tight=True)
#ax5.xaxis.set_ticklabels([])
#ax6.yaxis.set_label_position("left")
#ax6.yaxis.tick_left()
plt.ylabel('Gaussian width (m)', color='k', fontsize=12)
plt.xlabel('Along track distance (km)', color='k', fontsize=12)
ax5.yaxis.set_label_position("right")
ax5.yaxis.tick_right()
ax5.tick_params(axis='y', colors='k', labelsize=12)
plt.ylim(0,np.ceil(width_max))
#plt.xlim(ax5.get_xlim())
plt.xlim([sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1]])


plt.subplots_adjust(left=0.06, right=0.93, top=0.96, bottom=0.07, wspace=0, hspace=0)
### Save plot to file
plt.savefig(figPath+'/S-2_overlay_'+date+beam_name+relStr+str(Slim)+str(Nlim)+example+'radiometry.png', dpi=300)


