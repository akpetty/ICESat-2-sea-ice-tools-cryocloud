""" utils.py
    
    Random functions used by the various scripts and notebooks in this repo
    Initial code written by Alek Petty and Marco Bagnardi (02/01/2020) 

    Python dependencies:
        See below for the relevant module imports. More information on installation is given in the README file.

    Update history:
        02/01/2020: Version 1.
"""

import matplotlib, sys
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
import h5py
from glob import glob
#import seaborn as sns
from tqdm import tqdm
from astropy.time import Time
from datetime import datetime
import pyproj

def get_atl09_data_beam(ATL09, beam='gt1r', km_unit=False):
    """ Grab ATL09 data and put in pandas dataframe
    
    Args:
        ATL09 (hdf5 file): ATL09 file
        beamStr (str): ATLAS beam 
    
    Returns:
        ATL09_df: pandas dataframe

    Caveats:
        can loop through profiles if needed.
        ATL09 uses profiles (beam pair) not beam

    """

    if beam[2] == '1':
        profile = 'profile_1'
    elif beam[2] == '2':
        profile = 'profile_2'
    elif beam[2] == '3':
        profile = 'profile_3'
    
    # GPS time from start of ATLAS epoch
    delta_time_this=ATL09[profile + '/high_rate/delta_time'][:]

    # Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL09['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # add the time from GPS epoch to the ATLAS epoch
    gpsseconds = atlas_epoch + delta_time_this

    # Along-track distance
    dist_this = ATL09[profile + '/high_rate/prof_dist_x'][:]
    # Convert to km
    if km_unit:
        dist_this = dist_this * 0.001
    
    # Use astropy to convert GPS time to UTC datetime
    tiso=Time(gpsseconds,format='gps').utc.datetime

    # Beam dataframe
    ATL09_df = pd.DataFrame({
                                #
                                'lats': ATL09[profile + '/high_rate/latitude'][:],
                                'lons': ATL09[profile + '/high_rate/longitude'][:],
                                'delta_time': ATL09[profile + '/high_rate/delta_time'][:],
                                'along_dist': dist_this,
                                'layer_flag': ATL09[profile + '/high_rate/layer_flag'][:],
                                'cloud_flag_asr': ATL09[profile + '/high_rate/cloud_flag_asr'][:],
                                'cloud_flag_atm': ATL09[profile + '/high_rate/cloud_flag_atm'][:],
                                'utc_time':tiso,
                                })

    
    #Replace invalid value with NaN
    #ATL09_df = ATL09_df.replace(np.max(ATL09_df['LayerFlag']), np.nan)

    return ATL09_df

def get_atl07_data_beam_extra(ATL07, beamStr='', beamNum=0, km_unit=False, max_segment_length=400, min_segment_length=2):
    """ Grab ATL07 data and put in pandas dataframe, add some extra info for diagnostics
    
    Args:
        ATL07 (hdf5 file): ATL07 file
        beamStr (str or int): ATLAS ground track
        beamNum (int): ATLAS beam number
        km_unit (boleen): if true convert distance to km
    
    Returns:
        pd_data: pandas dataframe

    Caveats:
        need to add data gap condtion but not sure how..

    """

    # Find out the spacecraft orientation to know what beam string to choose
    orientation_flag=ATL07['orbit_info']['sc_orient'][0]
    print('orientation flag:', orientation_flag)
    if (orientation_flag==0):
        print('Backward orientation')
        beamStrs=['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
                
    elif (orientation_flag==1):
        print('Forward orientation')
        beamStrs=['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
        
    elif (orientation_flag==2):
        print('Transitioning, do not use for science!')
    
    if beamNum>0:
        # set ground track based on supplied beam number (if zero dont do this and use supplied ground track)
        beamStr=beamStrs[beamNum-1]
        print(beamStr)

    if (orientation_flag==0):
        #backward
        if (beamStr[-1]=='l'):
            beamStrength='strong'
        else:
            beamStrength='weak'
    elif (orientation_flag==1):
        #forward
        if (beamStr[-1]=='l'):
            beamStrength='weak'
        else:
            beamStrength='strong'
    else:
        beamStrength='NA'
        print('no clear spaceraft orientation, probably transitioning')


    # grab location data and along-track distance
    lon_this = ATL07[beamStr + '/sea_ice_segments/longitude'][:]
    lat_this = ATL07[beamStr + '/sea_ice_segments/latitude'][:]
    dist_this = ATL07[beamStr + '/sea_ice_segments/seg_dist_x'][:]

    # Convert to km
    if km_unit:
        dist_this = dist_this * 0.001

    # Height ID
    height_id_this = ATL07[beamStr + '/sea_ice_segments/height_segment_id'][:]

    # Height
    height_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_height'][:]

    # Height
    quality_flag_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_fit_quality_flag'][:]
    
    # Gaussian width
    height_width = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_w_gaussian'][:]
    print(height_width)
    # Set Gaussian width to nan where unphysically high
    height_width[np.where(height_width>10)]=np.nan
    # Seg Type
    seg_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_type'][:]
    # SSH_Flag
    ssh_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_ssh_flag'][:]
    # segment length
    seglength_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_length_seg'][:]
    # photon rate
    photon_rate = ATL07[beamStr + '/sea_ice_segments/stats/photon_rate'][:]
    
    # background rate
    bground_rate = ATL07[beamStr + '/sea_ice_segments/stats/backgr_r_200'][:]

    # background rate normalized
    bground_rate_norm = ATL07[beamStr + '/sea_ice_segments/stats/background_r_norm'][:]

    # Convert background rates to MHz
    bground_rate = bground_rate * 0.000001
    bground_rate_norm = bground_rate_norm* 0.000001

    # Solar elevation
    solar_elevation = ATL07[beamStr + '/sea_ice_segments/geolocation/solar_elevation'][:]

    # ice_conc
    #ice_conc = ATL07[beamStr + '/sea_ice_segments/stats/ice_conc'][:]
     # GPS time from start of ATLAS epoch
    delta_time_this=ATL07[beamStr + '/sea_ice_segments/delta_time'][:]

    # Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL07['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # add the time from GPS epoch to the ATLAS epoch
    gpsseconds = atlas_epoch + delta_time_this

    ## Use astropy to convert GPS time to UTC time
    tiso=Time(gpsseconds,format='gps').utc.datetime

     # seg stats
    seg_h_mean_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_mean_h'][:]
    seg_h_width_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_w'][:]
    seg_surf_err_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_surface_error_est'][:]

    #segs_used_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_n_pulse_seg_used'][:]
    #segs_total_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_n_pulse_seg'][:]

    # Cloud stats
    cloud_flag_asr_this = ATL07[beamStr + '/sea_ice_segments/stats/cloud_flag_asr'][:]
    cloud_flag_atm_this = ATL07[beamStr + '/sea_ice_segments/stats/cloud_flag_atm'][:]
    layer_flag_this = ATL07[beamStr + '/sea_ice_segments/stats/layer_flag'][:]

    # Beam gain
    #beam_gain = ATL07[beamStr + '/ancillary_data/beam_gain'][:]

    # Assign to pandas dataframe
    pd_data = pd.DataFrame({'lons': lon_this, 'lats': lat_this, 'along_dist':dist_this, 'gpsseconds':gpsseconds, 'utc_time':tiso, 
            'height':height_this, 'height_width':height_width, 'quality_flag':quality_flag_this, 'height_id':height_id_this, 'seg_type':seg_this, 'ssh_flag':ssh_this, 'seg_length': seglength_this, 
            'photon_rate':photon_rate, 'bground_rate':bground_rate, 'bground_rate_norm':bground_rate_norm, 'solar_elevation':solar_elevation,
            'cloud_flag_asr':cloud_flag_asr_this, 'cloud_flag_atm':cloud_flag_atm_this, 'layer_flag':layer_flag_this,
            'seg_h_mean':seg_h_mean_this, 'seg_h_width':seg_h_width_this, 'seg_surf_err':seg_surf_err_this})
            #'segs_used':segs_used_this, 'segs_total':segs_total_this,})
    

    print('seg limits:', np.amin(pd_data['seg_length'].values), np.amax(pd_data['seg_length'].values))
    print('lon limits:', np.amin(pd_data['lons'].values), np.amax(pd_data['lons'].values))
    print('lat limits:', np.amin(pd_data['lats'].values), np.amax(pd_data['lats'].values))

    # Filter data
    pd_data= pd_data.dropna() 
    pd_data = pd_data[pd_data['seg_length']<max_segment_length]
    pd_data = pd_data[pd_data['seg_length']>min_segment_length]
    pd_data = pd_data[pd_data['height']<1e6]
    pd_data = pd_data[pd_data['seg_h_mean']<1e6]


    print(ATL07, beamStr, beamStrength)

    return pd_data, beamStr, beamStrength

def get_atl10_freeboards_rel006(f1, beam, km_unit=False, max_segment_length=400, min_segment_length=2):
    """ Pandas/numpy ATL10 reader including changes for rel006
    Original function written by Alek Petty, June 2018 (alek.a.petty@nasa.gov)

    I've picked out the variables from ATL10 I think are of most interest to sea ice users, 
    but by no means is this an exhastive list. 
    See the xarray or dictionary readers to load in the more complete ATL10 dataset
    or explore the hdf5 files themselves (I like using the app Panpoly for this) to see what else you might want
    
    Args:
        fileT (str): File path of the ATL10 dataset
        beamStr (str): ICESat-2 beam (the number is the pair, r=strong, l=weak)
        maxFreeboard (float): maximum freeboard (meters)
    
    Optional args:
        epsg_string (str): EPSG string for projecting data (default of 3411 north polar stereo)

    returns:
        pandas dataframe
        
    Versions:
        v1: June 2018
        v2: June 2020 - cleaned things up, changed the time function slightly to be consistent with Ellen's ATL07 reader.

    """

    try:
        freeboard=f1[beam]['freeboard_segment']['beam_fb_height'][:]
    except:
        print ('no valid freeboard data for that beam')
        return None, beam

    # Freeboards
    #freeboard=f1[beam]['freeboard_beam_segment']['beam_freeboard']['beam_fb_height'][:]
    ssh_flag=f1[beam]['freeboard_segment']['heights']['height_segment_ssh_flag'][:]
    seg_type_flag=f1[beam]['freeboard_segment']['heights']['height_segment_type'][:]
    
    
    # Along track distance from start of granule
    dist = f1[beam]['freeboard_segment']['seg_dist_x'][:]
    # Convert to km
    if km_unit:
        dist = dist * 0.001
        
    # Height segment ID (10 km segments)
    height_segment_id=f1[beam]['freeboard_segment']['height_segment_id'][:]


    lead_loc=f1[beam +'/leads/ssh_ndx'][:]
    #lead_id = height_segment_id[lead_loc]
    lead_flag = np.zeros((np.size(height_segment_id)))
    # ADD SIZE??
    lead_flag[lead_loc]=1
    
    lons=f1[beam]['freeboard_segment']['longitude'][:]
    lats=f1[beam]['freeboard_segment']['latitude'][:]
    
    seg_length = f1[beam]['freeboard_segment']['heights']['height_segment_length_seg'][:]

    
     # GPS time from start of ATLAS epoch
    delta_time=f1[beam]['freeboard_segment']['delta_time'][:]

    # Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=f1['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # add the time from GPS epoch to the ATLAS epoch
    gps_seconds = atlas_epoch + delta_time

    ## Use astropy to convert GPS time to UTC time
    tiso=Time(gps_seconds,format='gps').utc.datetime
    
    pd_data = pd.DataFrame({'freeboard':freeboard, 'lons':lons, 'lats':lats, 'delta_time':delta_time,
                      'along_dist':dist, 'height_segment_id': height_segment_id, 'lead_flag':lead_flag, 
                        'ssh_flag':ssh_flag, 'seg_type_flag':seg_type_flag, 'seg_length': seg_length, 'utc_time': tiso})


    print('seg limits:', np.amin(pd_data['seg_length'].values), np.amax(pd_data['seg_length'].values))
    print('lon limits:', np.amin(pd_data['lons'].values), np.amax(pd_data['lons'].values))
    print('lat limits:', np.amin(pd_data['lats'].values), np.amax(pd_data['lats'].values))

    # Filter data
    pd_data= pd_data.dropna() 
    pd_data = pd_data[pd_data['seg_length']<max_segment_length]
    pd_data = pd_data[pd_data['seg_length']>min_segment_length]
    pd_data = pd_data[pd_data['freeboard']<1e6]

    
    return pd_data, beam
    
def get_atl10_data_beam_extra(ATL10, beamStr, km_unit=False, max_segment_length=400, min_segment_length=2):
    """ Grab ATL10 data and put in pandas dataframe, add some extra info for diagnostics
    
    Args:
        ATL10 (hdf5 file): ATL10 file
        beamStr (str or int): ATLAS ground track
        beamNum (int): ATLAS beam number
        km_unit (boleen): if true convert distance to km

    Returns:
        pd_data: pandas dataframe

    Caveats:
        need to add data gap condtion but not sure how..

    """

    ## Read ATL10 parameters of interest


    freeboard = ATL10[beamStr + '/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:]
    lons = ATL10[beamStr + '/freeboard_beam_segment/beam_freeboard/longitude'][:]
    lats = ATL10[beamStr + '/freeboard_beam_segment/beam_freeboard/latitude'][:]
    dist = ATL10[beamStr + '/freeboard_beam_segment/beam_freeboard/seg_dist_x'][:]
    delta_time=ATL10[beamStr + '/freeboard_beam_segment/beam_freeboard/delta_time'][:]
        

    # Convert to km
    if km_unit:
        dist = dist * 0.001

    height_segment_id = ATL10[beamStr + '/freeboard_beam_segment/beam_freeboard/height_segment_id'][:]
    ssh_flag=ATL10[beamStr + '/freeboard_beam_segment/height_segments/height_segment_ssh_flag'][:]
    seg_type_flag=ATL10[beamStr + '/freeboard_beam_segment/height_segments/height_segment_type'][:]
    seg_length=ATL10[beamStr + '/freeboard_beam_segment/height_segments/height_segment_length_seg'][:]
    
    lead_loc=ATL10[beamStr +'/leads/ssh_ndx'][:]
    #lead_id = height_segment_id[lead_loc]
    lead_flag = np.zeros((np.size(height_segment_id)))
    # ADD SIZE??
    lead_flag[lead_loc]=1

    # Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL10['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # add the time from GPS epoch to the ATLAS epoch
    gpsseconds = atlas_epoch + delta_time
    
    # Use astropy to convert GPS time to UTC datetime
    tiso=Time(gpsseconds,format='gps').utc.datetime


    pd_data = pd.DataFrame({'freeboard':freeboard, 'lons':lons, 'lats':lats, 'delta_time':delta_time,
                      'along_dist':dist, 'height_segment_id': height_segment_id, 'lead_flag':lead_flag, 
                        'ssh_flag':ssh_flag, 'seg_type_flag':seg_type_flag, 'seg_length': seg_length, 'utc_time': tiso})


    print('seg limits:', np.amin(pd_data['seg_length'].values), np.amax(pd_data['seg_length'].values))
    print('lon limits:', np.amin(pd_data['lons'].values), np.amax(pd_data['lons'].values))
    print('lat limits:', np.amin(pd_data['lats'].values), np.amax(pd_data['lats'].values))

    # Filter data
    pd_data= pd_data.dropna() 
    pd_data = pd_data[pd_data['seg_length']<max_segment_length]
    pd_data = pd_data[pd_data['seg_length']>min_segment_length]
    pd_data = pd_data[pd_data['freeboard']<1e6]

    return pd_data, beamStr

def get_cpom_chords(data_path, hem, start_date='201811', end_date='201904', version='v2', min=0.5, max=50):
    """ Get CPOM chord length data """
    lats=np.array([])
    lons=np.array([])
    chords=np.array([])
    if hem=='nh':
        files_all=glob(data_path+version+'/*'+'arco_size.dat')
    elif hem=='sh':
        files_all=glob(data_path+version+'/*'+'anto_size.dat')
    date_range = pd.date_range(start_date,end_date, freq='M').strftime("%Y%m").tolist()
    
    cpom_files_range = [s for s in files_all for sub in date_range if sub in s ]
    print(cpom_files_range)

    for file in cpom_files_range:
        data = np.loadtxt(file)
        lats=np.hstack([lats, data[:, 0]])
        lons=np.hstack([lons,data[:, 1]])
        chords=np.hstack([chords, data[:, 2]])
    
    mask=np.where((chords<max) & (chords>min) & ~np.isnan(chords))
    chords=chords[mask]
    lats=lats[mask]
    lons=lons[mask]
    
    return lats, lons, chords

def calc_leadheight_thresh(data, idxT, height_var='height', seg_h_width_var='seg_h_width', seg_h_mean_var='seg_h_mean', seg_surf_err_var='seg_surf_err',percentile_thresh=20):
    """ Calculate the local lead height threshold
    See Section 4.1.1.17 of the Sea Ice ATBD
    
    Args:
        data (pandas dataframe): original ATL07 data
        idxT (list): section index
        percentile_thresh (int): snow density variable of choosing
    
    Returns:
        lead_htmax: lead height threshold

    Caveats:
         - I'm not sure this is exactly the same as the ATBD, still checking this. 
         - need to add a data gap condtion but not sure how
         - not sure what to do if there are no smooth heights
         - not sure exactly why the histogram mean is used here
         - should we take the max error in the section or use the height segment error for each possible lead segment? i.e should this return one number for the section or an array. 

    """

    # 1. Calculate the nth percentile (percentile_thresh) of the heights in this section
    
    height_lower_pct = np.percentile(data.iloc[idxT][height_var].values, percentile_thresh)
    #print(height_thresh)
    
    # 2. Calculate the minimum height from the smooth height histogram means and the expected height error

    # Find the indices of the 'smooth' or 'good' height retrievals
    idx_emin=np.where(data.iloc[idxT][seg_h_width_var].values<0.13)
    # Find the minimum histogram height
    #print('num smooth segs', np.size(idx_emin))
    if len(idx_emin)>2:
        emin=np.nanmin(data.iloc[idxT][seg_h_mean_var].values[idx_emin])
    else:
        # Set to a very low value, basically meaning we pick the sorted height percentile instead
        emin=-999
    
    # Add the expected error to this 
    min_height_plus_error = emin + (2.*data.iloc[idxT][seg_surf_err_var]) 
    min_height_plus_error=min_height_plus_error[np.isfinite(min_height_plus_error)]

    # Find the max value (I'm not sure this is the same as the ATBD)
    min_height_plus_error_max=np.nanmax(min_height_plus_error)
    lead_htmax = max(min_height_plus_error_max, height_lower_pct)
    
    #print('lead htmax:', lead_htmax, min_height_plus_error_max, height_lower_pct, emin)
    return lead_htmax


def calc_chords(leads, pd_data, version_str, minsize=3, minlength=60., maxdiff=300, max_chord_length=50000, footprint_size=11):
    """ Finding chord lengths
    
     Args:
        leads (numpy array): array of along track distance with lead segments set to negative numbers
        pd_data (dataframe): pandas dataframe of ATL07 data
        minsize (int): minimum number of segments for a valid chord
        minlength (int): minimum length for a valid chord
        maxdiff (int): maximum segment gap within a chord
        max_chord_length (int): maximum chord length (set to nan if higher than this)
        footprint_size (int): footprint size, add this to the chord length estimate (half each side)

    To do:
        Since switching to pandas it makes sense to use groupby to group the floes instead of this ad-hoc attempt.
        Replace loops with 
    
    Version 2 (03/10/21): 
        Switched to just not processing a chord instead of setting to nan if greater than 50km. 
        added segment ids (start and end) to match against thickness data.    

    """

    splitdata=np.split(leads, np.where(leads< -100)[0])

    # Should really just split/group the entire pandas dataframe!
    splitdata_lon=np.split(pd_data.lons.values, np.where(leads< -100)[0])
    splitdata_lat=np.split(pd_data.lats.values, np.where(leads< -100)[0])
    splitdata_time=np.split(pd_data.gpsseconds.values, np.where(leads< -100)[0])
    splitdata_id=np.split(pd_data.height_id.values, np.where(leads< -100)[0])
    splitdata_atrack=np.split(pd_data.along_dist.values, np.where(leads< -100)[0])


    lengths=[]
    groups=[]
    lons=[]
    atrack=[]
    lats=[]
    time=[]
    start_id=[]
    end_id=[]
    num_segs=[]

    #print('first group from granule', splitdata[0])
    print('processing chords...')
    for x in np.arange(np.size(splitdata)):
        # start at index 1 as index 0 is the lead
        # If it's a chord (i.e. positive values)
        if np.mean(splitdata[x][1:])>-900:
            # If greater that minimum number of segments
            if np.size(splitdata[x][1:]) >= minsize:
                #If greater than minimum length 
                 if (splitdata[x][-1]-splitdata[x][1]) >= minlength:
                    #If less than minimum length 
                    if (splitdata[x][-1]-splitdata[x][1]) <= max_chord_length:
                        #If no segment gap greater than max gap 
                        if (np.max(np.diff(splitdata[x][1:]))<maxdiff):
                            
                            #print('height ID:',splitdata_id[x])
                            #print('height ID s', splitdata_id[x][1])
                            #print('height ID e', splitdata_id[x][-1])
                            #print('group atrk:',splitdata[x])


                            groups.append(splitdata[x][1:])
                            lengths.append(splitdata[x][-1]-splitdata[x][1]+footprint_size)
                            lons.append(np.mean(splitdata_lon[x][1:]))
                            lats.append(np.mean(splitdata_lat[x][1:]))
                            time.append(np.mean(splitdata_time[x][1:]))
                            atrack.append(np.mean(splitdata_atrack[x][1:]))
                            
                            start_id.append(splitdata_id[x][1])
                            end_id.append(splitdata_id[x][-1])
                            num_segs.append(np.size(splitdata[x][1:]))
                        else:
                            print('dropped floe group as max segment separation = ', np.max(np.diff(splitdata[x][1:])))
                    else:
                        print('dropped floe group as floe length = ', splitdata[x][-1]-splitdata[x][1])    


    print('Num chords:', np.size(lengths))
    print('Mean lengths:', np.mean(lengths))
    print('Mean start_id:', np.mean(start_id))
    print('Mean end id:', np.mean(end_id))
    print('Mean num segs:', np.mean(num_segs))

    pd_data_chords = pd.DataFrame({'chord_lengths': lengths, 'lons': lons, 'lats':lats, 'along_track_distance':atrack, 'time': time, 'start_chord_id': start_id, 'end_chord_id': end_id, 'num_segs':num_segs})
    #print(pd_data_chords)

    # nan chord lengths > 50 km
    #pd_data_chords.loc[pd_data_chords['lengths']>max_chord_length, ['lengths', 'num_segs', 'start_chord_id']] = np.nan

    return pd_data_chords, groups


def calc_chords_old(leads, pd_data, version_str, minsize=3, minlength=60., maxdiff=300, max_chord_length=50000, footprint_size=11):
    """ Finding chord lengths
    
     Args:
        leads (numpy array): array of along track distance with lead segments set to negative numbers
        pd_data (dataframe): pandas dataframe of ATL07 data
        minsize (int): minimum number of segments for a valid chord
        minlength (int): minimum length for a valid chord
        maxdiff (int): maximum segment gap within a chord
        max_chord_length (int): maximum chord length (set to nan if higher than this)
        footprint_size (int): footprint size, add this to the chord length estimate (half each side)

    To do:
        Since switching to pandas it makes sense to use groupby to group the floes instead of this ad-hoc attempt.
    
    """

    splitdata=np.split(leads, np.where(leads< -100)[0])

    splitdata_lon=np.split(pd_data.lons, np.where(leads< -100)[0])
    splitdata_lat=np.split(pd_data.lats, np.where(leads< -100)[0])
    splitdata_time=np.split(pd_data.gpsseconds, np.where(leads< -100)[0])


    lengths=[]
    groups=[]
    lons=[]
    lats=[]
    time=[]

    #print('first group from granule', splitdata[0])
    print('processing chords...')
    for x in np.arange(np.size(splitdata)):
        # start at index 1 as index 0 is the lead
        # If it's a chord (i.e. positive values)
        if np.mean(splitdata[x][1:])>-900:
            # If greater that minimum number of segments
            if np.size(splitdata[x][1:]) >= minsize:
                #If greater than minimum length 
                 if (splitdata[x][-1]-splitdata[x][1]) >= minlength:
                    #If no segment gap greater than max gap 
                    if (np.max(np.diff(splitdata[x][1:]))<maxdiff):
                        
                        groups.append(splitdata[x][1:])
                        lengths.append(splitdata[x][-1]-splitdata[x][1]+footprint_size)
                        lons.append(np.mean(splitdata_lon[x][1:]))
                        lats.append(np.mean(splitdata_lat[x][1:]))
                        time.append(np.mean(splitdata_time[x][1:]))
                    else:
                        print('dropped floe group as max segment separation = ', np.max(np.diff(splitdata[x][1:])))
                        
    #print('last group from granule', splitdata[x])

    pd_data_chords = pd.DataFrame({'lengths': lengths, 'lons': lons, 'lats':lats, 'time': time})
    #print(pd_data_chords)

    # nan chord lenghts > 50 km
    pd_data_chords.loc[pd_data_chords['lengths']>max_chord_length, ['lengths']] = np.nan

    return pd_data_chords, groups



def get_chords_2versions(data, km_unit=True, lon_var = 'lons', lat_var='lats', alongtrack_var='along_dist', sshvar = 'ssh_flag', segtype_var = 'seg_type', height_var='height', percentileThresh=20, percentile_section_size=10000):
    """ Calculate two versions of chord lengths from ATL07

    Use the SSH flag and also a variable height filter to seperate and calculate chord lengths
    
    Args:
        data (pandas dataframe): data
        percentileThresh (int): local perncetile height threshold
        percentile_section_size (int): size of sections for calculating local height threshold (in meters)
    
    Returns:
        version 1 and version 2 chord length data
    
    Caveats:
        should probably split this into two separate functions (v1 and v2)
    """

    if km_unit:
        data=data.copy()
        # convert along track to meters
        data[alongtrack_var] = data[alongtrack_var].multiply(1000.)

    
    lead_indexv1=np.copy(data[alongtrack_var]) 
    lead_indexv2=np.copy(data[alongtrack_var]) 
    #print('along min/max:', np.amin(lead_indexv2), np.amax(lead_indexv2))

    lead_indexv1[np.where(data[sshvar]>0.5)]=-999
    
    # loop over 10 km sections to do the percentile height /seg type classification as in ATL10 processing
    for x in range(int(data[alongtrack_var].iloc[0]), np.int(data[alongtrack_var].iloc[-1]), percentile_section_size):
        #print(x, percentile_section_size)
        #print('x:', x)
        idxT=np.where((lead_indexv2>x)&(lead_indexv2<(x+percentile_section_size)))[0]
        #print('idx:', idxT)
        if (np.size(idxT)>10):
            # require at least 10 segments

            height_thresh = calc_leadheight_thresh(data, idxT, percentile_thresh=percentileThresh)
            #print(np.percentile(data[height_var].iloc[idxT], percentileThresh))
            idxM=np.where((data[segtype_var].iloc[idxT]>1.5 ) & (data[segtype_var].iloc[idxT]<5.5 ) & (data[height_var].iloc[idxT]<height_thresh))
            #print(idxT[0][0]+idxM)
            lead_indexv2[idxT[0]+idxM]=-999

        else:
            continue
    
    # Find raw floe lengths
    pddata_v1, groupsv1 = calc_chords(lead_indexv1, data, 'v1')
    pddata_v2, groupsv2 = calc_chords(lead_indexv2, data, 'v2')
            
    return pddata_v1, pddata_v2, lead_indexv1, lead_indexv2, groupsv1, groupsv2

def get_chords_3versions(data, km_unit=True, lon_var = 'lons', lat_var='lats', alongtrack_var='along_dist', sshvar = 'ssh_flag', segtype_var = 'seg_type', height_var='height', percentileThresh=20, percentile_section_size=10000):
    """ Calculate three versions of chord lengths from ATL07

    Use the SSH flag and also a variable height filter to seperate and calculate chord lengths
    
    Args:
        data (pandas dataframe): data
        percentileThresh (int): local perncetile height threshold
        percentile_section_size (int): size of sections for calculating local height threshold (in meters)
    
    Returns:
        version 1, 2 and 3 chord length data
    
    Caveats:
        should probably split this into two separate functions (v1 and v2)
    """

    if km_unit:
        data=data.copy()
        # convert along track to meters
        data[alongtrack_var] = data[alongtrack_var].multiply(1000.)

    
    lead_indexv1=np.copy(data[alongtrack_var]) 
    lead_indexv2=np.copy(data[alongtrack_var]) 
    lead_indexv3=np.copy(data[alongtrack_var]) 
    #print('along min/max:', np.amin(lead_indexv2), np.amax(lead_indexv2))

    lead_indexv1[np.where(data[sshvar]>0.5)]=-999
    
    # loop over 10 km sections to do the percentile height /seg type classification as in ATL10 processing
    for x in range(int(data[alongtrack_var].iloc[0]), np.int(data[alongtrack_var].iloc[-1]), percentile_section_size):

        idxT=np.where((lead_indexv2>x)&(lead_indexv2<(x+percentile_section_size)))[0]
        
        if (np.size(idxT)>10):
            # require at least 10 segments

            height_thresh = calc_leadheight_thresh(data, idxT, percentile_thresh=percentileThresh)
            
            # Include only specular leads
            idxM_v2=np.where((data[segtype_var].iloc[idxT]>1.5 ) & (data[segtype_var].iloc[idxT]<5.5 ) & (data[height_var].iloc[idxT]<height_thresh))
            
            # Include dark leads too
            idxM_v3=np.where((data[segtype_var].iloc[idxT]>1.5 ) & (data[segtype_var].iloc[idxT]<9.5 ) & (data[height_var].iloc[idxT]<height_thresh))
            #print(idxT[0][0]+idxM)
            lead_indexv2[idxT[0]+idxM_v2]=-999
            lead_indexv3[idxT[0]+idxM_v3]=-999

        else:
            continue
    
    # Find raw floe lengths
    pddata_v1, groupsv1 = calc_chords(lead_indexv1, data, 'v1')
    pddata_v2, groupsv2 = calc_chords(lead_indexv2, data, 'v2')
    pddata_v3, groupsv3 = calc_chords(lead_indexv3, data, 'v3')
            
    return pddata_v1, pddata_v2, pddata_v3, lead_indexv1, lead_indexv2, lead_indexv3, groupsv1, groupsv2, groupsv3


def reset_matplotlib():
    """ Reset matplotlib to a common default. """
    
    # Force agg backend.
    plt.switch_backend('agg')
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    plt.rcParams['ytick.major.size'] = 2
    plt.rcParams['axes.linewidth'] = .6
    plt.rcParams['lines.linewidth'] = .6
    plt.rcParams['patch.linewidth'] = .6
    plt.rcParams['ytick.labelsize']=8
    plt.rcParams['xtick.labelsize']=8
    plt.rcParams['legend.fontsize']=9
    plt.rcParams['font.size']=9
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def get_beam_config(ATL07file):
    """ Get the beam string configuration for the spacecraft orientation of the given file"""
    
    orientation_flag=ATL07file['orbit_info']['sc_orient'][0]
    print('orientation flag:', orientation_flag)
    
    if (orientation_flag==0):
        print('Backward orientation')
        beamStrs=['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
                
    elif (orientation_flag==1):
        print('Forward orientation')
        beamStrs=['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
        
    elif (orientation_flag==2):
        print('Transitioning, do not use for science!')
    
    return beamStrs

def create_grid_nsidc_arctic(epsg_string='3411', nx=304, ny=448, leftx=-3837500, dxRes=25000, uppery=5837500, dyRes=25000):
    """ Use pyproj to generate a grid covering the given domain (defined by the projection and the corner lat/lons)"""

    crs = pyproj.CRS.from_string("epsg:"+epsg_string)
    p=pyproj.Proj(crs)
    
    print(dxRes, dyRes)

    x=leftx+dxRes*np.indices((ny,nx),np.float32)[1]
    y=uppery-dxRes*np.indices((ny,nx),np.float32)[0]
    lons, lats = p(x, y, inverse=True)
    
    return x, y, lats, lons, p


def create_grid_nsidc_antarctic(epsg_string='3412', nx=316, ny=332, leftx=-3937500, dxRes=25000, uppery=4337500, dyRes=25000):
    """ Use pyproj to generate a grid covering the given domain (defined by the projection and the corner lat/lons)"""

    crs = pyproj.CRS.from_string("epsg:"+epsg_string)
    p=pyproj.Proj(crs)
    
    print(dxRes, dyRes)

    x=leftx+dxRes*np.indices((ny,nx),np.float32)[1]
    y=uppery-dxRes*np.indices((ny,nx),np.float32)[0]
    lons, lats = p(x, y, inverse=True)
    
    return x, y, lats, lons, p


def bin_data_nsidc(xpts, ypts, var, xptsG, yptsG, dx):
    """ Bin data using numpy histogram 2d
    
    Adapted for the NSIDC grid which has orgin in the top left corner.

    """
    xptsG2=np.flipud(xptsG)
    yptsG2=np.flipud(yptsG)

    xbins=xptsG2[0]-(dx/2)
    ybins=yptsG2[:, 0]-(dx/2)
    # add on one more bin in each direction
    xbins2=np.append(xbins, xbins[-1]+dx)
    ybins2=np.append(ybins, ybins[-1]+dx)
    print('binning..')
    print(xbins2.shape)
    print(ybins2.shape)
    z, _, _ = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2), weights=var)
    counts, _, _ = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2))

    varG = z / counts
    
    # Need to flip the arrayback as we did this to make the y axis go from negative to positive, then we need to transpose because of weirdness with histogram2d
    varG=np.flipud(varG.T)
    counts=np.flipud(counts.T)

    return varG, counts

def get_region_mask(ancDataPath, mplot, xypts_return=0):
    header = 300
    datatype='uint8'
    file_mask = ancDataPath+'/region_n.msk'
    
    #1 Non-regional ocean  
    #2 Sea of Okhotsk 
    #3 Bering Sea  
    #4 Hudson Bay 
    #5 Baffin Bay/Davis Strait/Labrador Sea    
    #6 Greenland Sea   Bellingshausen 
    #7 Kara and Barents Seas

    #8 - Arctic Ocean
    #9 - Canadian Archipelago
    #10 - Gulf of St Lawrence
    #11 - Land

    fd = open(file_mask, 'rb')
    region_mask = np.fromfile(file=fd, dtype=datatype)
    region_mask = np.reshape(region_mask[header:], [448, 304])

    if (xypts_return==1):
        mask_latf = open(ancDataPath+'/psn25lats_v3.dat', 'rb')
        mask_lonf = open(ancDataPath+'/psn25lons_v3.dat', 'rb')
        lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
        lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

        xpts, ypts = mplot(lons_mask, lats_mask)

        return region_mask, xpts, ypts
    else:
        return region_mask


def get_atl07_data_beam(ATL07, beamStr):
    """ Gran ATL07 data and put in pandas dataframer
    
    Args:
        ATL07 (hdf5 file): ATL07 file
        beamStr (str): ATLAS beam 
    
    Returns:
        pd_data: pandas dataframe

    Caveats:
        need to add data gap condtion but not sure how..

    """

    # grab location data and along-track distance
    lon_this = ATL07[beamStr + '/sea_ice_segments/longitude'][:]
    lat_this = ATL07[beamStr + '/sea_ice_segments/latitude'][:]
    dist_this = ATL07[beamStr + '/sea_ice_segments/seg_dist_x'][:]
    # Height
    height_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_height'][:]
    # Seg Type
    seg_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_type'][:]
    # SSH_Flag
    ssh_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_ssh_flag'][:]
    # segment length
    seglength_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_length_seg'][:]
    # ice_conc
    ice_conc = ATL07[beamStr + '/sea_ice_segments/stats/ice_conc'][:]
    
    # GPS time from start of ATLAS epoch
    delta_time_this=ATL07[beamStr + '/sea_ice_segments/delta_time'][:]

    # Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL07['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # add the time from GPS epoch to the ATLAS epoch
    gpsseconds = atlas_epoch + delta_time_this

    ## Use astropy to convert GPS time to UTC time
    #tiso=Time(gps_seconds,format='gps').utc.datetime

     # seg stats
    seg_h_mean_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_mean_h'][:]
    seg_h_width_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_w'][:]
    seg_surf_err_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_surface_error_est'][:]

    # Assign to pandas dataframe
    pd_data = pd.DataFrame({'lons': lon_this, 'lats': lat_this, 'along_dist':dist_this, 'seg_length': seglength_this, 
            'seg_type':seg_this, 'height':height_this, 'ssh_flag':ssh_this, 'seg_h_mean':seg_h_mean_this, 
            'seg_h_width':seg_h_width_this, 'seg_surf_err':seg_surf_err_this, 'gpsseconds':gpsseconds, 'ice_conc':ice_conc})
    

    print('seg limits:', np.amin(pd_data['seg_length'].values), np.amax(pd_data['seg_length'].values))
    print('lon limits:', np.amin(pd_data['lons'].values), np.amax(pd_data['lons'].values))
    print('lat limits:', np.amin(pd_data['lats'].values), np.amax(pd_data['lats'].values))

    # Filter data
    pd_data= pd_data.dropna() 
    pd_data = pd_data[pd_data['seg_length']<200]
    pd_data = pd_data[pd_data['seg_length']>2]
    pd_data = pd_data[pd_data['height']<1e6]
    pd_data = pd_data[pd_data['seg_h_mean']<1e6]

    return pd_data


def get_region_mask_s(ancDataPath, mplot, xypts_return=0):
    header = 300
    datatype='uint8'
    file_mask = ancDataPath+'/region_s.msk'
    

    #2 Weddell Sea
    #3 Indian Ocean  
    #4 Pacific Ocean
    #5 Ross Sea
    #6 Bellingshausen Amundsen Sea

    #11 - Land
    #12 - Coast

    fd = open(file_mask, 'rb')
    region_mask = np.fromfile(file=fd, dtype=datatype)
    region_mask = np.reshape(region_mask[header:], [332, 316])

    if (xypts_return==1):
        mask_latf = open(ancDataPath+'/pss25lats_v3.dat', 'rb')
        mask_lonf = open(ancDataPath+'/pss25lons_v3.dat', 'rb')
        lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [332, 316])
        lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [332, 316])

        xpts, ypts = mplot(lons_mask, lats_mask)

        return region_mask, xpts, ypts
    else:
        return region_mask
