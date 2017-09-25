#%% Python file set-up
# 
import os
# NetCDF data handler
import netCDF4 as nc
# numerical computing package
import numpy as np
# basemap toolkit to plot maps
from mpl_toolkits.basemap import Basemap
# command style functions that make matplotlib work like MATLAB
import matplotlib.pyplot as plt
# matplotlib
import matplotlib
# find nearest value
def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()
#
import sys

os.chdir('/Users/earl/PROJECT2')
figures_path = '/Users/earl/Dropbox/PROJECT2/figures/'
data_path = '/Users/earl/Dropbox/Data/PROJECT2/'

Mon = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', \
       'Oct', 'Nov', 'Dec']
MMM = ['JFM', 'AMJ', 'JAS', 'OND']

scriptname = os.path.basename(sys.argv[0])[:-3]

print('OK, all set')


#%% Load data tau_x CT
# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_0' + '.nc', 'r')

# same for temperature (1,2,3)
temp_JFM_CT_0 = nc_fid.variables['temp'][0,:,:,:]

# get dimensions
lat = nc_fid.variables['yt_ocean'][:]
lon1 = nc_fid.variables['xt_ocean'][:]
depth_0 = nc_fid.variables['st_ocean'][:]

#% Feb
nc_fid = nc.Dataset(data_path + 'KDS75/temp_AMJ_0' + '.nc', 'r')
temp_AMJ_CT_0 = nc_fid.variables['temp'][0,:,:,:]
# Mar
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JAS_0' + '.nc', 'r')
temp_JAS_CT_0 = nc_fid.variables['temp'][0,:,:,:]
# MMM average
nc_fid = nc.Dataset(data_path + 'KDS75/temp_OND_0' + '.nc', 'r')
temp_OND_CT_0 = nc_fid.variables['temp'][0,:,:,:]

print('CT temp 0 OK')


#%% temp 1
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_1' + '.nc', 'r')
temp_JFM_CT_1 = nc_fid.variables['temp'][0,:,:,:]
depth_1 = nc_fid.variables['st_ocean'][:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_AMJ_1' + '.nc', 'r')
temp_AMJ_CT_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JAS_1' + '.nc', 'r')
temp_JAS_CT_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_OND_1' + '.nc', 'r')
temp_OND_CT_1 = nc_fid.variables['temp'][0,:,:,:]
print('CT temp 1 OK')


#%% temp 2
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_2' + '.nc', 'r')
temp_JFM_CT_2 = nc_fid.variables['temp'][0,:,:,:]
depth_2 = nc_fid.variables['st_ocean'][:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_AMJ_2' + '.nc', 'r')
temp_AMJ_CT_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JAS_2' + '.nc', 'r')
temp_JAS_CT_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_OND_2' + '.nc', 'r')
temp_OND_CT_2 = nc_fid.variables['temp'][0,:,:,:]
print('CT temp 2 OK')


#%% temp 3
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_3' + '.nc', 'r')
temp_JFM_CT_3 = nc_fid.variables['temp'][0,:,:,:]
depth_3 = nc_fid.variables['st_ocean'][:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_AMJ_3' + '.nc', 'r')
temp_AMJ_CT_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JAS_3' + '.nc', 'r')
temp_JAS_CT_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/temp_OND_3' + '.nc', 'r')
temp_OND_CT_3 = nc_fid.variables['temp'][0,:,:,:]
print('CT temp 3 OK')


#%% temp UP
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JFM_0' + '.nc', 'r')
temp_JFM_UP_0 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_AMJ_0' + '.nc', 'r')
temp_AMJ_UP_0 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JAS_0' + '.nc', 'r')
temp_JAS_UP_0 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_OND_0' + '.nc', 'r')
temp_OND_UP_0 = nc_fid.variables['temp'][0,:,:,:]
print('UP temp 0 OK')

#% temp 1
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JFM_1' + '.nc', 'r')
temp_JFM_UP_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_AMJ_1' + '.nc', 'r')
temp_AMJ_UP_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JAS_1' + '.nc', 'r')
temp_JAS_UP_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_OND_1' + '.nc', 'r')
temp_OND_UP_1 = nc_fid.variables['temp'][0,:,:,:]
print('UP temp 1 OK')

#% temp 2
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JFM_2' + '.nc', 'r')
temp_JFM_UP_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_AMJ_2' + '.nc', 'r')
temp_AMJ_UP_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JAS_2' + '.nc', 'r')
temp_JAS_UP_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_OND_2' + '.nc', 'r')
temp_OND_UP_2 = nc_fid.variables['temp'][0,:,:,:]
print('UP temp 2 OK')

#% temp 3
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JFM_3' + '.nc', 'r')
temp_JFM_UP_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_AMJ_3' + '.nc', 'r')
temp_AMJ_UP_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_JAS_3' + '.nc', 'r')
temp_JAS_UP_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/temp_OND_3' + '.nc', 'r')
temp_OND_UP_3 = nc_fid.variables['temp'][0,:,:,:]
print('UP temp 3 OK')



#%% temp PI
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JFM_0' + '.nc', 'r')
temp_JFM_PI_0 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_AMJ_0' + '.nc', 'r')
temp_AMJ_PI_0 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JAS_0' + '.nc', 'r')
temp_JAS_PI_0 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_OND_0' + '.nc', 'r')
temp_OND_PI_0 = nc_fid.variables['temp'][0,:,:,:]
print('PI temp 0 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JFM_1' + '.nc', 'r')
temp_JFM_PI_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_AMJ_1' + '.nc', 'r')
temp_AMJ_PI_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JAS_1' + '.nc', 'r')
temp_JAS_PI_1 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_OND_1' + '.nc', 'r')
temp_OND_PI_1 = nc_fid.variables['temp'][0,:,:,:]
print('PI temp 1 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JFM_2' + '.nc', 'r')
temp_JFM_PI_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_AMJ_2' + '.nc', 'r')
temp_AMJ_PI_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JAS_2' + '.nc', 'r')
temp_JAS_PI_2 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_OND_2' + '.nc', 'r')
temp_OND_PI_2 = nc_fid.variables['temp'][0,:,:,:]
print('PI temp 2 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JFM_3' + '.nc', 'r')
temp_JFM_PI_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_AMJ_3' + '.nc', 'r')
temp_AMJ_PI_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_JAS_3' + '.nc', 'r')
temp_JAS_PI_3 = nc_fid.variables['temp'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/temp_OND_3' + '.nc', 'r')
temp_OND_PI_3 = nc_fid.variables['temp'][0,:,:,:]
print('PI temp 3 OK')



#%% salt CT UP PI
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JFM_0' + '.nc', 'r')
salt_JFM_CT_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_AMJ_0' + '.nc', 'r')
salt_AMJ_CT_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JAS_0' + '.nc', 'r')
salt_JAS_CT_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_OND_0' + '.nc', 'r')
salt_OND_CT_0 = nc_fid.variables['salt'][0,:,:,:]
print('CT salt 1 OK')
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JFM_1' + '.nc', 'r')
salt_JFM_CT_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_AMJ_1' + '.nc', 'r')
salt_AMJ_CT_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JAS_1' + '.nc', 'r')
salt_JAS_CT_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_OND_1' + '.nc', 'r')
salt_OND_CT_1 = nc_fid.variables['salt'][0,:,:,:]
print('CT salt 1 OK')
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JFM_2' + '.nc', 'r')
salt_JFM_CT_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_AMJ_2' + '.nc', 'r')
salt_AMJ_CT_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JAS_2' + '.nc', 'r')
salt_JAS_CT_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_OND_2' + '.nc', 'r')
salt_OND_CT_2 = nc_fid.variables['salt'][0,:,:,:]
print('CT salt 2 OK')
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JFM_3' + '.nc', 'r')
salt_JFM_CT_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_AMJ_3' + '.nc', 'r')
salt_AMJ_CT_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_JAS_3' + '.nc', 'r')
salt_JAS_CT_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/salt_OND_3' + '.nc', 'r')
salt_OND_CT_3 = nc_fid.variables['salt'][0,:,:,:]
print('CT salt 3 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JFM_0' + '.nc', 'r')
salt_JFM_UP_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_AMJ_0' + '.nc', 'r')
salt_AMJ_UP_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JAS_0' + '.nc', 'r')
salt_JAS_UP_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_OND_0' + '.nc', 'r')
salt_OND_UP_0 = nc_fid.variables['salt'][0,:,:,:]
print('UP salt 0 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JFM_1' + '.nc', 'r')
salt_JFM_UP_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_AMJ_1' + '.nc', 'r')
salt_AMJ_UP_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JAS_1' + '.nc', 'r')
salt_JAS_UP_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_OND_1' + '.nc', 'r')
salt_OND_UP_1 = nc_fid.variables['salt'][0,:,:,:]
print('UP salt 1 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JFM_2' + '.nc', 'r')
salt_JFM_UP_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_AMJ_2' + '.nc', 'r')
salt_AMJ_UP_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JAS_2' + '.nc', 'r')
salt_JAS_UP_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_OND_2' + '.nc', 'r')
salt_OND_UP_2 = nc_fid.variables['salt'][0,:,:,:]
print('UP salt 2 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JFM_3' + '.nc', 'r')
salt_JFM_UP_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_AMJ_3' + '.nc', 'r')
salt_AMJ_UP_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_JAS_3' + '.nc', 'r')
salt_JAS_UP_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/salt_OND_3' + '.nc', 'r')
salt_OND_UP_3 = nc_fid.variables['salt'][0,:,:,:]
print('UP salt 3 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JFM_0' + '.nc', 'r')
salt_JFM_PI_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_AMJ_0' + '.nc', 'r')
salt_AMJ_PI_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JAS_0' + '.nc', 'r')
salt_JAS_PI_0 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_OND_0' + '.nc', 'r')
salt_OND_PI_0 = nc_fid.variables['salt'][0,:,:,:]
print('PI salt 0 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JFM_1' + '.nc', 'r')
salt_JFM_PI_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_AMJ_1' + '.nc', 'r')
salt_AMJ_PI_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JAS_1' + '.nc', 'r')
salt_JAS_PI_1 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_OND_1' + '.nc', 'r')
salt_OND_PI_1 = nc_fid.variables['salt'][0,:,:,:]
print('PI salt 1 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JFM_2' + '.nc', 'r')
salt_JFM_PI_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_AMJ_2' + '.nc', 'r')
salt_AMJ_PI_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JAS_2' + '.nc', 'r')
salt_JAS_PI_2 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_OND_2' + '.nc', 'r')
salt_OND_PI_2 = nc_fid.variables['salt'][0,:,:,:]
print('PI salt 2 OK')
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JFM_3' + '.nc', 'r')
salt_JFM_PI_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_AMJ_3' + '.nc', 'r')
salt_AMJ_PI_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_JAS_3' + '.nc', 'r')
salt_JAS_PI_3 = nc_fid.variables['salt'][0,:,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/salt_OND_3' + '.nc', 'r')
salt_OND_PI_3 = nc_fid.variables['salt'][0,:,:,:]
print('PI salt 3 OK')



#%% re-organise longitudes
lon_less_m70 = lon1 < -70
lon_less_m70_flipped = np.flipud(lon_less_m70)

temp_JFM_fl_CT_0 = np.ma.empty((15,988,3600))
temp_JFM_fl_CT_0[:,:,lon_less_m70_flipped] = temp_JFM_CT_0[:,:,lon_less_m70];
temp_JFM_fl_CT_0[:,:,~lon_less_m70_flipped] = temp_JFM_CT_0[:,:,~lon_less_m70];
temp_AMJ_fl_CT_0 = np.ma.empty((15,988,3600))
temp_AMJ_fl_CT_0[:,:,lon_less_m70_flipped] = temp_AMJ_CT_0[:,:,lon_less_m70];
temp_AMJ_fl_CT_0[:,:,~lon_less_m70_flipped] = temp_AMJ_CT_0[:,:,~lon_less_m70];
temp_JAS_fl_CT_0 = np.ma.empty((15,988,3600))
temp_JAS_fl_CT_0[:,:,lon_less_m70_flipped] = temp_JAS_CT_0[:,:,lon_less_m70];
temp_JAS_fl_CT_0[:,:,~lon_less_m70_flipped] = temp_JAS_CT_0[:,:,~lon_less_m70];
temp_OND_fl_CT_0 = np.ma.empty((15,988,3600))
temp_OND_fl_CT_0[:,:,lon_less_m70_flipped] = temp_OND_CT_0[:,:,lon_less_m70];
temp_OND_fl_CT_0[:,:,~lon_less_m70_flipped] = temp_OND_CT_0[:,:,~lon_less_m70];

#%%
temp_JFM_fl_UP_0 = np.ma.empty((15,988,3600))
temp_JFM_fl_UP_0[:,:,lon_less_m70_flipped] = temp_JFM_UP_0[:,:,lon_less_m70];
temp_JFM_fl_UP_0[:,:,~lon_less_m70_flipped] = temp_JFM_UP_0[:,:,~lon_less_m70];
temp_AMJ_fl_UP_0 = np.ma.empty((15,988,3600))
temp_AMJ_fl_UP_0[:,:,lon_less_m70_flipped] = temp_AMJ_UP_0[:,:,lon_less_m70];
temp_AMJ_fl_UP_0[:,:,~lon_less_m70_flipped] = temp_AMJ_UP_0[:,:,~lon_less_m70];
temp_JAS_fl_UP_0 = np.ma.empty((15,988,3600))
temp_JAS_fl_UP_0[:,:,lon_less_m70_flipped] = temp_JAS_UP_0[:,:,lon_less_m70];
temp_JAS_fl_UP_0[:,:,~lon_less_m70_flipped] = temp_JAS_UP_0[:,:,~lon_less_m70];
temp_OND_fl_UP_0 = np.ma.empty((15,988,3600))
temp_OND_fl_UP_0[:,:,lon_less_m70_flipped] = temp_OND_UP_0[:,:,lon_less_m70];
temp_OND_fl_UP_0[:,:,~lon_less_m70_flipped] = temp_OND_UP_0[:,:,~lon_less_m70];

#%%
temp_JFM_fl_PI_0 = np.ma.empty((15,988,3600))
temp_JFM_fl_PI_0[:,:,lon_less_m70_flipped] = temp_JFM_PI_0[:,:,lon_less_m70];
temp_JFM_fl_PI_0[:,:,~lon_less_m70_flipped] = temp_JFM_PI_0[:,:,~lon_less_m70];
temp_AMJ_fl_PI_0 = np.ma.empty((15,988,3600))
temp_AMJ_fl_PI_0[:,:,lon_less_m70_flipped] = temp_AMJ_PI_0[:,:,lon_less_m70];
temp_AMJ_fl_PI_0[:,:,~lon_less_m70_flipped] = temp_AMJ_PI_0[:,:,~lon_less_m70];
temp_JAS_fl_PI_0 = np.ma.empty((15,988,3600))
temp_JAS_fl_PI_0[:,:,lon_less_m70_flipped] = temp_JAS_PI_0[:,:,lon_less_m70];
temp_JAS_fl_PI_0[:,:,~lon_less_m70_flipped] = temp_JAS_PI_0[:,:,~lon_less_m70];
temp_OND_fl_PI_0 = np.ma.empty((15,988,3600))
temp_OND_fl_PI_0[:,:,lon_less_m70_flipped] = temp_OND_PI_0[:,:,lon_less_m70];
temp_OND_fl_PI_0[:,:,~lon_less_m70_flipped] = temp_OND_PI_0[:,:,~lon_less_m70];


#%% 1
temp_JFM_fl_CT_1 = np.ma.empty((15,988,3600))
temp_JFM_fl_CT_1[:,:,lon_less_m70_flipped] = temp_JFM_CT_1[:,:,lon_less_m70];
temp_JFM_fl_CT_1[:,:,~lon_less_m70_flipped] = temp_JFM_CT_1[:,:,~lon_less_m70];
temp_AMJ_fl_CT_1 = np.ma.empty((15,988,3600))
temp_AMJ_fl_CT_1[:,:,lon_less_m70_flipped] = temp_AMJ_CT_1[:,:,lon_less_m70];
temp_AMJ_fl_CT_1[:,:,~lon_less_m70_flipped] = temp_AMJ_CT_1[:,:,~lon_less_m70];
temp_JAS_fl_CT_1 = np.ma.empty((15,988,3600))
temp_JAS_fl_CT_1[:,:,lon_less_m70_flipped] = temp_JAS_CT_1[:,:,lon_less_m70];
temp_JAS_fl_CT_1[:,:,~lon_less_m70_flipped] = temp_JAS_CT_1[:,:,~lon_less_m70];
temp_OND_fl_CT_1 = np.ma.empty((15,988,3600))
temp_OND_fl_CT_1[:,:,lon_less_m70_flipped] = temp_OND_CT_1[:,:,lon_less_m70];
temp_OND_fl_CT_1[:,:,~lon_less_m70_flipped] = temp_OND_CT_1[:,:,~lon_less_m70];
temp_JFM_fl_UP_1 = np.ma.empty((15,988,3600))
temp_JFM_fl_UP_1[:,:,lon_less_m70_flipped] = temp_JFM_UP_1[:,:,lon_less_m70];
temp_JFM_fl_UP_1[:,:,~lon_less_m70_flipped] = temp_JFM_UP_1[:,:,~lon_less_m70];
temp_AMJ_fl_UP_1 = np.ma.empty((15,988,3600))
temp_AMJ_fl_UP_1[:,:,lon_less_m70_flipped] = temp_AMJ_UP_1[:,:,lon_less_m70];
temp_AMJ_fl_UP_1[:,:,~lon_less_m70_flipped] = temp_AMJ_UP_1[:,:,~lon_less_m70];
temp_JAS_fl_UP_1 = np.ma.empty((15,988,3600))
temp_JAS_fl_UP_1[:,:,lon_less_m70_flipped] = temp_JAS_UP_1[:,:,lon_less_m70];
temp_JAS_fl_UP_1[:,:,~lon_less_m70_flipped] = temp_JAS_UP_1[:,:,~lon_less_m70];
temp_OND_fl_UP_1 = np.ma.empty((15,988,3600))
temp_OND_fl_UP_1[:,:,lon_less_m70_flipped] = temp_OND_UP_1[:,:,lon_less_m70];
temp_OND_fl_UP_1[:,:,~lon_less_m70_flipped] = temp_OND_UP_1[:,:,~lon_less_m70];
temp_JFM_fl_PI_1 = np.ma.empty((15,988,3600))
temp_JFM_fl_PI_1[:,:,lon_less_m70_flipped] = temp_JFM_PI_1[:,:,lon_less_m70];
temp_JFM_fl_PI_1[:,:,~lon_less_m70_flipped] = temp_JFM_PI_1[:,:,~lon_less_m70];
temp_AMJ_fl_PI_1 = np.ma.empty((15,988,3600))
temp_AMJ_fl_PI_1[:,:,lon_less_m70_flipped] = temp_AMJ_PI_1[:,:,lon_less_m70];
temp_AMJ_fl_PI_1[:,:,~lon_less_m70_flipped] = temp_AMJ_PI_1[:,:,~lon_less_m70];
temp_JAS_fl_PI_1 = np.ma.empty((15,988,3600))
temp_JAS_fl_PI_1[:,:,lon_less_m70_flipped] = temp_JAS_PI_1[:,:,lon_less_m70];
temp_JAS_fl_PI_1[:,:,~lon_less_m70_flipped] = temp_JAS_PI_1[:,:,~lon_less_m70];
temp_OND_fl_PI_1 = np.ma.empty((15,988,3600))
temp_OND_fl_PI_1[:,:,lon_less_m70_flipped] = temp_OND_PI_1[:,:,lon_less_m70];
temp_OND_fl_PI_1[:,:,~lon_less_m70_flipped] = temp_OND_PI_1[:,:,~lon_less_m70];


#%% 2
temp_JFM_fl_CT_2 = np.ma.empty((15,988,3600))
temp_JFM_fl_CT_2[:,:,lon_less_m70_flipped] = temp_JFM_CT_2[:,:,lon_less_m70];
temp_JFM_fl_CT_2[:,:,~lon_less_m70_flipped] = temp_JFM_CT_2[:,:,~lon_less_m70];
temp_AMJ_fl_CT_2 = np.ma.empty((15,988,3600))
temp_AMJ_fl_CT_2[:,:,lon_less_m70_flipped] = temp_AMJ_CT_2[:,:,lon_less_m70];
temp_AMJ_fl_CT_2[:,:,~lon_less_m70_flipped] = temp_AMJ_CT_2[:,:,~lon_less_m70];
temp_JAS_fl_CT_2 = np.ma.empty((15,988,3600))
temp_JAS_fl_CT_2[:,:,lon_less_m70_flipped] = temp_JAS_CT_2[:,:,lon_less_m70];
temp_JAS_fl_CT_2[:,:,~lon_less_m70_flipped] = temp_JAS_CT_2[:,:,~lon_less_m70];
temp_OND_fl_CT_2 = np.ma.empty((15,988,3600))
temp_OND_fl_CT_2[:,:,lon_less_m70_flipped] = temp_OND_CT_2[:,:,lon_less_m70];
temp_OND_fl_CT_2[:,:,~lon_less_m70_flipped] = temp_OND_CT_2[:,:,~lon_less_m70];
temp_JFM_fl_UP_2 = np.ma.empty((15,988,3600))
temp_JFM_fl_UP_2[:,:,lon_less_m70_flipped] = temp_JFM_UP_2[:,:,lon_less_m70];
temp_JFM_fl_UP_2[:,:,~lon_less_m70_flipped] = temp_JFM_UP_2[:,:,~lon_less_m70];
temp_AMJ_fl_UP_2 = np.ma.empty((15,988,3600))
temp_AMJ_fl_UP_2[:,:,lon_less_m70_flipped] = temp_AMJ_UP_2[:,:,lon_less_m70];
temp_AMJ_fl_UP_2[:,:,~lon_less_m70_flipped] = temp_AMJ_UP_2[:,:,~lon_less_m70];
temp_JAS_fl_UP_2 = np.ma.empty((15,988,3600))
temp_JAS_fl_UP_2[:,:,lon_less_m70_flipped] = temp_JAS_UP_2[:,:,lon_less_m70];
temp_JAS_fl_UP_2[:,:,~lon_less_m70_flipped] = temp_JAS_UP_2[:,:,~lon_less_m70];
temp_OND_fl_UP_2 = np.ma.empty((15,988,3600))
temp_OND_fl_UP_2[:,:,lon_less_m70_flipped] = temp_OND_UP_2[:,:,lon_less_m70];
temp_OND_fl_UP_2[:,:,~lon_less_m70_flipped] = temp_OND_UP_2[:,:,~lon_less_m70];
temp_JFM_fl_PI_2 = np.ma.empty((15,988,3600))
temp_JFM_fl_PI_2[:,:,lon_less_m70_flipped] = temp_JFM_PI_2[:,:,lon_less_m70];
temp_JFM_fl_PI_2[:,:,~lon_less_m70_flipped] = temp_JFM_PI_2[:,:,~lon_less_m70];
temp_AMJ_fl_PI_2 = np.ma.empty((15,988,3600))
temp_AMJ_fl_PI_2[:,:,lon_less_m70_flipped] = temp_AMJ_PI_2[:,:,lon_less_m70];
temp_AMJ_fl_PI_2[:,:,~lon_less_m70_flipped] = temp_AMJ_PI_2[:,:,~lon_less_m70];
temp_JAS_fl_PI_2 = np.ma.empty((15,988,3600))
temp_JAS_fl_PI_2[:,:,lon_less_m70_flipped] = temp_JAS_PI_2[:,:,lon_less_m70];
temp_JAS_fl_PI_2[:,:,~lon_less_m70_flipped] = temp_JAS_PI_2[:,:,~lon_less_m70];
temp_OND_fl_PI_2 = np.ma.empty((15,988,3600))
temp_OND_fl_PI_2[:,:,lon_less_m70_flipped] = temp_OND_PI_2[:,:,lon_less_m70];
temp_OND_fl_PI_2[:,:,~lon_less_m70_flipped] = temp_OND_PI_2[:,:,~lon_less_m70];


#%% 3
temp_JFM_fl_CT_3 = np.ma.empty((15,988,3600))
temp_JFM_fl_CT_3[:,:,lon_less_m70_flipped] = temp_JFM_CT_3[:,:,lon_less_m70];
temp_JFM_fl_CT_3[:,:,~lon_less_m70_flipped] = temp_JFM_CT_3[:,:,~lon_less_m70];
temp_AMJ_fl_CT_3 = np.ma.empty((15,988,3600))
temp_AMJ_fl_CT_3[:,:,lon_less_m70_flipped] = temp_AMJ_CT_3[:,:,lon_less_m70];
temp_AMJ_fl_CT_3[:,:,~lon_less_m70_flipped] = temp_AMJ_CT_3[:,:,~lon_less_m70];
temp_JAS_fl_CT_3 = np.ma.empty((15,988,3600))
temp_JAS_fl_CT_3[:,:,lon_less_m70_flipped] = temp_JAS_CT_3[:,:,lon_less_m70];
temp_JAS_fl_CT_3[:,:,~lon_less_m70_flipped] = temp_JAS_CT_3[:,:,~lon_less_m70];
temp_OND_fl_CT_3 = np.ma.empty((15,988,3600))
temp_OND_fl_CT_3[:,:,lon_less_m70_flipped] = temp_OND_CT_3[:,:,lon_less_m70];
temp_OND_fl_CT_3[:,:,~lon_less_m70_flipped] = temp_OND_CT_3[:,:,~lon_less_m70];
temp_JFM_fl_UP_3 = np.ma.empty((15,988,3600))
temp_JFM_fl_UP_3[:,:,lon_less_m70_flipped] = temp_JFM_UP_3[:,:,lon_less_m70];
temp_JFM_fl_UP_3[:,:,~lon_less_m70_flipped] = temp_JFM_UP_3[:,:,~lon_less_m70];
temp_AMJ_fl_UP_3 = np.ma.empty((15,988,3600))
temp_AMJ_fl_UP_3[:,:,lon_less_m70_flipped] = temp_AMJ_UP_3[:,:,lon_less_m70];
temp_AMJ_fl_UP_3[:,:,~lon_less_m70_flipped] = temp_AMJ_UP_3[:,:,~lon_less_m70];
temp_JAS_fl_UP_3 = np.ma.empty((15,988,3600))
temp_JAS_fl_UP_3[:,:,lon_less_m70_flipped] = temp_JAS_UP_3[:,:,lon_less_m70];
temp_JAS_fl_UP_3[:,:,~lon_less_m70_flipped] = temp_JAS_UP_3[:,:,~lon_less_m70];
temp_OND_fl_UP_3 = np.ma.empty((15,988,3600))
temp_OND_fl_UP_3[:,:,lon_less_m70_flipped] = temp_OND_UP_3[:,:,lon_less_m70];
temp_OND_fl_UP_3[:,:,~lon_less_m70_flipped] = temp_OND_UP_3[:,:,~lon_less_m70];
temp_JFM_fl_PI_3 = np.ma.empty((15,988,3600))
temp_JFM_fl_PI_3[:,:,lon_less_m70_flipped] = temp_JFM_PI_3[:,:,lon_less_m70];
temp_JFM_fl_PI_3[:,:,~lon_less_m70_flipped] = temp_JFM_PI_3[:,:,~lon_less_m70];
temp_AMJ_fl_PI_3 = np.ma.empty((15,988,3600))
temp_AMJ_fl_PI_3[:,:,lon_less_m70_flipped] = temp_AMJ_PI_3[:,:,lon_less_m70];
temp_AMJ_fl_PI_3[:,:,~lon_less_m70_flipped] = temp_AMJ_PI_3[:,:,~lon_less_m70];
temp_JAS_fl_PI_3 = np.ma.empty((15,988,3600))
temp_JAS_fl_PI_3[:,:,lon_less_m70_flipped] = temp_JAS_PI_3[:,:,lon_less_m70];
temp_JAS_fl_PI_3[:,:,~lon_less_m70_flipped] = temp_JAS_PI_3[:,:,~lon_less_m70];
temp_OND_fl_PI_3 = np.ma.empty((15,988,3600))
temp_OND_fl_PI_3[:,:,lon_less_m70_flipped] = temp_OND_PI_3[:,:,lon_less_m70];
temp_OND_fl_PI_3[:,:,~lon_less_m70_flipped] = temp_OND_PI_3[:,:,~lon_less_m70];


#%% salt
salt_JFM_fl_CT_0 = np.ma.empty((15,988,3600))
salt_JFM_fl_CT_0[:,:,lon_less_m70_flipped] = salt_JFM_CT_0[:,:,lon_less_m70];
salt_JFM_fl_CT_0[:,:,~lon_less_m70_flipped] = salt_JFM_CT_0[:,:,~lon_less_m70];
salt_AMJ_fl_CT_0 = np.ma.empty((15,988,3600))
salt_AMJ_fl_CT_0[:,:,lon_less_m70_flipped] = salt_AMJ_CT_0[:,:,lon_less_m70];
salt_AMJ_fl_CT_0[:,:,~lon_less_m70_flipped] = salt_AMJ_CT_0[:,:,~lon_less_m70];
salt_JAS_fl_CT_0 = np.ma.empty((15,988,3600))
salt_JAS_fl_CT_0[:,:,lon_less_m70_flipped] = salt_JAS_CT_0[:,:,lon_less_m70];
salt_JAS_fl_CT_0[:,:,~lon_less_m70_flipped] = salt_JAS_CT_0[:,:,~lon_less_m70];
salt_OND_fl_CT_0 = np.ma.empty((15,988,3600))
salt_OND_fl_CT_0[:,:,lon_less_m70_flipped] = salt_OND_CT_0[:,:,lon_less_m70];
salt_OND_fl_CT_0[:,:,~lon_less_m70_flipped] = salt_OND_CT_0[:,:,~lon_less_m70];
salt_JFM_fl_UP_0 = np.ma.empty((15,988,3600))
salt_JFM_fl_UP_0[:,:,lon_less_m70_flipped] = salt_JFM_UP_0[:,:,lon_less_m70];
salt_JFM_fl_UP_0[:,:,~lon_less_m70_flipped] = salt_JFM_UP_0[:,:,~lon_less_m70];
salt_AMJ_fl_UP_0 = np.ma.empty((15,988,3600))
salt_AMJ_fl_UP_0[:,:,lon_less_m70_flipped] = salt_AMJ_UP_0[:,:,lon_less_m70];
salt_AMJ_fl_UP_0[:,:,~lon_less_m70_flipped] = salt_AMJ_UP_0[:,:,~lon_less_m70];
salt_JAS_fl_UP_0 = np.ma.empty((15,988,3600))
salt_JAS_fl_UP_0[:,:,lon_less_m70_flipped] = salt_JAS_UP_0[:,:,lon_less_m70];
salt_JAS_fl_UP_0[:,:,~lon_less_m70_flipped] = salt_JAS_UP_0[:,:,~lon_less_m70];
salt_OND_fl_UP_0 = np.ma.empty((15,988,3600))
salt_OND_fl_UP_0[:,:,lon_less_m70_flipped] = salt_OND_UP_0[:,:,lon_less_m70];
salt_OND_fl_UP_0[:,:,~lon_less_m70_flipped] = salt_OND_UP_0[:,:,~lon_less_m70];
salt_JFM_fl_PI_0 = np.ma.empty((15,988,3600))
salt_JFM_fl_PI_0[:,:,lon_less_m70_flipped] = salt_JFM_PI_0[:,:,lon_less_m70];
salt_JFM_fl_PI_0[:,:,~lon_less_m70_flipped] = salt_JFM_PI_0[:,:,~lon_less_m70];
salt_AMJ_fl_PI_0 = np.ma.empty((15,988,3600))
salt_AMJ_fl_PI_0[:,:,lon_less_m70_flipped] = salt_AMJ_PI_0[:,:,lon_less_m70];
salt_AMJ_fl_PI_0[:,:,~lon_less_m70_flipped] = salt_AMJ_PI_0[:,:,~lon_less_m70];
salt_JAS_fl_PI_0 = np.ma.empty((15,988,3600))
salt_JAS_fl_PI_0[:,:,lon_less_m70_flipped] = salt_JAS_PI_0[:,:,lon_less_m70];
salt_JAS_fl_PI_0[:,:,~lon_less_m70_flipped] = salt_JAS_PI_0[:,:,~lon_less_m70];
salt_OND_fl_PI_0 = np.ma.empty((15,988,3600))
salt_OND_fl_PI_0[:,:,lon_less_m70_flipped] = salt_OND_PI_0[:,:,lon_less_m70];
salt_OND_fl_PI_0[:,:,~lon_less_m70_flipped] = salt_OND_PI_0[:,:,~lon_less_m70];
salt_JFM_fl_CT_1 = np.ma.empty((15,988,3600))
salt_JFM_fl_CT_1[:,:,lon_less_m70_flipped] = salt_JFM_CT_1[:,:,lon_less_m70];
salt_JFM_fl_CT_1[:,:,~lon_less_m70_flipped] = salt_JFM_CT_1[:,:,~lon_less_m70];
salt_AMJ_fl_CT_1 = np.ma.empty((15,988,3600))
salt_AMJ_fl_CT_1[:,:,lon_less_m70_flipped] = salt_AMJ_CT_1[:,:,lon_less_m70];
salt_AMJ_fl_CT_1[:,:,~lon_less_m70_flipped] = salt_AMJ_CT_1[:,:,~lon_less_m70];
salt_JAS_fl_CT_1 = np.ma.empty((15,988,3600))
salt_JAS_fl_CT_1[:,:,lon_less_m70_flipped] = salt_JAS_CT_1[:,:,lon_less_m70];
salt_JAS_fl_CT_1[:,:,~lon_less_m70_flipped] = salt_JAS_CT_1[:,:,~lon_less_m70];
salt_OND_fl_CT_1 = np.ma.empty((15,988,3600))
salt_OND_fl_CT_1[:,:,lon_less_m70_flipped] = salt_OND_CT_1[:,:,lon_less_m70];
salt_OND_fl_CT_1[:,:,~lon_less_m70_flipped] = salt_OND_CT_1[:,:,~lon_less_m70];
salt_JFM_fl_UP_1 = np.ma.empty((15,988,3600))
salt_JFM_fl_UP_1[:,:,lon_less_m70_flipped] = salt_JFM_UP_1[:,:,lon_less_m70];
salt_JFM_fl_UP_1[:,:,~lon_less_m70_flipped] = salt_JFM_UP_1[:,:,~lon_less_m70];
salt_AMJ_fl_UP_1 = np.ma.empty((15,988,3600))
salt_AMJ_fl_UP_1[:,:,lon_less_m70_flipped] = salt_AMJ_UP_1[:,:,lon_less_m70];
salt_AMJ_fl_UP_1[:,:,~lon_less_m70_flipped] = salt_AMJ_UP_1[:,:,~lon_less_m70];
salt_JAS_fl_UP_1 = np.ma.empty((15,988,3600))
salt_JAS_fl_UP_1[:,:,lon_less_m70_flipped] = salt_JAS_UP_1[:,:,lon_less_m70];
salt_JAS_fl_UP_1[:,:,~lon_less_m70_flipped] = salt_JAS_UP_1[:,:,~lon_less_m70];
salt_OND_fl_UP_1 = np.ma.empty((15,988,3600))
salt_OND_fl_UP_1[:,:,lon_less_m70_flipped] = salt_OND_UP_1[:,:,lon_less_m70];
salt_OND_fl_UP_1[:,:,~lon_less_m70_flipped] = salt_OND_UP_1[:,:,~lon_less_m70];
salt_JFM_fl_PI_1 = np.ma.empty((15,988,3600))
salt_JFM_fl_PI_1[:,:,lon_less_m70_flipped] = salt_JFM_PI_1[:,:,lon_less_m70];
salt_JFM_fl_PI_1[:,:,~lon_less_m70_flipped] = salt_JFM_PI_1[:,:,~lon_less_m70];
salt_AMJ_fl_PI_1 = np.ma.empty((15,988,3600))
salt_AMJ_fl_PI_1[:,:,lon_less_m70_flipped] = salt_AMJ_PI_1[:,:,lon_less_m70];
salt_AMJ_fl_PI_1[:,:,~lon_less_m70_flipped] = salt_AMJ_PI_1[:,:,~lon_less_m70];
salt_JAS_fl_PI_1 = np.ma.empty((15,988,3600))
salt_JAS_fl_PI_1[:,:,lon_less_m70_flipped] = salt_JAS_PI_1[:,:,lon_less_m70];
salt_JAS_fl_PI_1[:,:,~lon_less_m70_flipped] = salt_JAS_PI_1[:,:,~lon_less_m70];
salt_OND_fl_PI_1 = np.ma.empty((15,988,3600))
salt_OND_fl_PI_1[:,:,lon_less_m70_flipped] = salt_OND_PI_1[:,:,lon_less_m70];
salt_OND_fl_PI_1[:,:,~lon_less_m70_flipped] = salt_OND_PI_1[:,:,~lon_less_m70];
salt_JFM_fl_CT_2 = np.ma.empty((15,988,3600))
salt_JFM_fl_CT_2[:,:,lon_less_m70_flipped] = salt_JFM_CT_2[:,:,lon_less_m70];
salt_JFM_fl_CT_2[:,:,~lon_less_m70_flipped] = salt_JFM_CT_2[:,:,~lon_less_m70];
salt_AMJ_fl_CT_2 = np.ma.empty((15,988,3600))
salt_AMJ_fl_CT_2[:,:,lon_less_m70_flipped] = salt_AMJ_CT_2[:,:,lon_less_m70];
salt_AMJ_fl_CT_2[:,:,~lon_less_m70_flipped] = salt_AMJ_CT_2[:,:,~lon_less_m70];
salt_JAS_fl_CT_2 = np.ma.empty((15,988,3600))
salt_JAS_fl_CT_2[:,:,lon_less_m70_flipped] = salt_JAS_CT_2[:,:,lon_less_m70];
salt_JAS_fl_CT_2[:,:,~lon_less_m70_flipped] = salt_JAS_CT_2[:,:,~lon_less_m70];
salt_OND_fl_CT_2 = np.ma.empty((15,988,3600))
salt_OND_fl_CT_2[:,:,lon_less_m70_flipped] = salt_OND_CT_2[:,:,lon_less_m70];
salt_OND_fl_CT_2[:,:,~lon_less_m70_flipped] = salt_OND_CT_2[:,:,~lon_less_m70];
salt_JFM_fl_UP_2 = np.ma.empty((15,988,3600))
salt_JFM_fl_UP_2[:,:,lon_less_m70_flipped] = salt_JFM_UP_2[:,:,lon_less_m70];
salt_JFM_fl_UP_2[:,:,~lon_less_m70_flipped] = salt_JFM_UP_2[:,:,~lon_less_m70];
salt_AMJ_fl_UP_2 = np.ma.empty((15,988,3600))
salt_AMJ_fl_UP_2[:,:,lon_less_m70_flipped] = salt_AMJ_UP_2[:,:,lon_less_m70];
salt_AMJ_fl_UP_2[:,:,~lon_less_m70_flipped] = salt_AMJ_UP_2[:,:,~lon_less_m70];
salt_JAS_fl_UP_2 = np.ma.empty((15,988,3600))
salt_JAS_fl_UP_2[:,:,lon_less_m70_flipped] = salt_JAS_UP_2[:,:,lon_less_m70];
salt_JAS_fl_UP_2[:,:,~lon_less_m70_flipped] = salt_JAS_UP_2[:,:,~lon_less_m70];
salt_OND_fl_UP_2 = np.ma.empty((15,988,3600))
salt_OND_fl_UP_2[:,:,lon_less_m70_flipped] = salt_OND_UP_2[:,:,lon_less_m70];
salt_OND_fl_UP_2[:,:,~lon_less_m70_flipped] = salt_OND_UP_2[:,:,~lon_less_m70];
salt_JFM_fl_PI_2 = np.ma.empty((15,988,3600))
salt_JFM_fl_PI_2[:,:,lon_less_m70_flipped] = salt_JFM_PI_2[:,:,lon_less_m70];
salt_JFM_fl_PI_2[:,:,~lon_less_m70_flipped] = salt_JFM_PI_2[:,:,~lon_less_m70];
salt_AMJ_fl_PI_2 = np.ma.empty((15,988,3600))
salt_AMJ_fl_PI_2[:,:,lon_less_m70_flipped] = salt_AMJ_PI_2[:,:,lon_less_m70];
salt_AMJ_fl_PI_2[:,:,~lon_less_m70_flipped] = salt_AMJ_PI_2[:,:,~lon_less_m70];
salt_JAS_fl_PI_2 = np.ma.empty((15,988,3600))
salt_JAS_fl_PI_2[:,:,lon_less_m70_flipped] = salt_JAS_PI_2[:,:,lon_less_m70];
salt_JAS_fl_PI_2[:,:,~lon_less_m70_flipped] = salt_JAS_PI_2[:,:,~lon_less_m70];
salt_OND_fl_PI_2 = np.ma.empty((15,988,3600))
salt_OND_fl_PI_2[:,:,lon_less_m70_flipped] = salt_OND_PI_2[:,:,lon_less_m70];
salt_OND_fl_PI_2[:,:,~lon_less_m70_flipped] = salt_OND_PI_2[:,:,~lon_less_m70];
salt_JFM_fl_CT_3 = np.ma.empty((15,988,3600))
salt_JFM_fl_CT_3[:,:,lon_less_m70_flipped] = salt_JFM_CT_3[:,:,lon_less_m70];
salt_JFM_fl_CT_3[:,:,~lon_less_m70_flipped] = salt_JFM_CT_3[:,:,~lon_less_m70];
salt_AMJ_fl_CT_3 = np.ma.empty((15,988,3600))
salt_AMJ_fl_CT_3[:,:,lon_less_m70_flipped] = salt_AMJ_CT_3[:,:,lon_less_m70];
salt_AMJ_fl_CT_3[:,:,~lon_less_m70_flipped] = salt_AMJ_CT_3[:,:,~lon_less_m70];
salt_JAS_fl_CT_3 = np.ma.empty((15,988,3600))
salt_JAS_fl_CT_3[:,:,lon_less_m70_flipped] = salt_JAS_CT_3[:,:,lon_less_m70];
salt_JAS_fl_CT_3[:,:,~lon_less_m70_flipped] = salt_JAS_CT_3[:,:,~lon_less_m70];
salt_OND_fl_CT_3 = np.ma.empty((15,988,3600))
salt_OND_fl_CT_3[:,:,lon_less_m70_flipped] = salt_OND_CT_3[:,:,lon_less_m70];
salt_OND_fl_CT_3[:,:,~lon_less_m70_flipped] = salt_OND_CT_3[:,:,~lon_less_m70];
salt_JFM_fl_UP_3 = np.ma.empty((15,988,3600))
salt_JFM_fl_UP_3[:,:,lon_less_m70_flipped] = salt_JFM_UP_3[:,:,lon_less_m70];
salt_JFM_fl_UP_3[:,:,~lon_less_m70_flipped] = salt_JFM_UP_3[:,:,~lon_less_m70];
salt_AMJ_fl_UP_3 = np.ma.empty((15,988,3600))
salt_AMJ_fl_UP_3[:,:,lon_less_m70_flipped] = salt_AMJ_UP_3[:,:,lon_less_m70];
salt_AMJ_fl_UP_3[:,:,~lon_less_m70_flipped] = salt_AMJ_UP_3[:,:,~lon_less_m70];
salt_JAS_fl_UP_3 = np.ma.empty((15,988,3600))
salt_JAS_fl_UP_3[:,:,lon_less_m70_flipped] = salt_JAS_UP_3[:,:,lon_less_m70];
salt_JAS_fl_UP_3[:,:,~lon_less_m70_flipped] = salt_JAS_UP_3[:,:,~lon_less_m70];
salt_OND_fl_UP_3 = np.ma.empty((15,988,3600))
salt_OND_fl_UP_3[:,:,lon_less_m70_flipped] = salt_OND_UP_3[:,:,lon_less_m70];
salt_OND_fl_UP_3[:,:,~lon_less_m70_flipped] = salt_OND_UP_3[:,:,~lon_less_m70];
salt_JFM_fl_PI_3 = np.ma.empty((15,988,3600))
salt_JFM_fl_PI_3[:,:,lon_less_m70_flipped] = salt_JFM_PI_3[:,:,lon_less_m70];
salt_JFM_fl_PI_3[:,:,~lon_less_m70_flipped] = salt_JFM_PI_3[:,:,~lon_less_m70];
salt_AMJ_fl_PI_3 = np.ma.empty((15,988,3600))
salt_AMJ_fl_PI_3[:,:,lon_less_m70_flipped] = salt_AMJ_PI_3[:,:,lon_less_m70];
salt_AMJ_fl_PI_3[:,:,~lon_less_m70_flipped] = salt_AMJ_PI_3[:,:,~lon_less_m70];
salt_JAS_fl_PI_3 = np.ma.empty((15,988,3600))
salt_JAS_fl_PI_3[:,:,lon_less_m70_flipped] = salt_JAS_PI_3[:,:,lon_less_m70];
salt_JAS_fl_PI_3[:,:,~lon_less_m70_flipped] = salt_JAS_PI_3[:,:,~lon_less_m70];
salt_OND_fl_PI_3 = np.ma.empty((15,988,3600))
salt_OND_fl_PI_3[:,:,lon_less_m70_flipped] = salt_OND_PI_3[:,:,lon_less_m70];
salt_OND_fl_PI_3[:,:,~lon_less_m70_flipped] = salt_OND_PI_3[:,:,~lon_less_m70];


#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
m = Basemap(projection='gall', llcrnrlat=-70,urcrnrlat=0.01,\
llcrnrlon=-70,urcrnrlon=290, resolution='h')
m.fix_aspect = True
m_aspect = m.aspect

lon = lon2 + 210


#%% plot surface maps: TEMPERATURE
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 6}) 

fig.set_size_inches(12, 5) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(431)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.03
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')

# meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
x, y = m(lons, lats)

# for quiver. data points
yy = np.arange(0, y.shape[0], 75)
xx = np.arange(0, x.shape[1], 125)
quiv_scale = 5

#
points = np.meshgrid(yy, xx)

# draw land outlines
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')

# filled contour plot
contf_lvls = np.arange(-2,31,1) # levels to show on colourbar
contf = m.contourf(x, y, temp_JFM_fl_CT_0[0,:,:], contf_lvls, cmap=cmap, extend='both')

# title ...
title_name = 'KDS75 temp JFM'
ax.set_title(title_name)

def round_to_base(x, base=5):
    return int(base * round(float(x) / base))
# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,0])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-80, 20, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])

#plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:2] \
#            + '_fig1_' + '.png', bbox_inches='tight', dpi=300)


#% 3000
ax = fig.add_subplot(434)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_AMJ_fl_CT_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 temp AMJ'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,1,0,0])
parallels = map(round_to_base, np.arange(-80, 20, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])


#% 2000
ax = fig.add_subplot(437)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_JAS_fl_CT_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 temp JAS'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[1,0,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])


# 1000
ax = fig.add_subplot(4,3,10)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_OND_fl_CT_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 temp OND'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])


#% 3000
ax = fig.add_subplot(432)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_JFM_fl_UP_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 UP temp JFM'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,1,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,0])

#% 3000
ax = fig.add_subplot(435)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_AMJ_fl_UP_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 UP temp AMJ'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,1,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,0])

# 2000
ax = fig.add_subplot(438)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_JAS_fl_UP_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 UP temp JAS'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[1,0,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,1])    


# 1000
ax = fig.add_subplot(4,3,11)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_OND_fl_UP_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 UP temp OND'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,1])


#% 3000
ax = fig.add_subplot(433)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_JFM_fl_PI_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 PI temp JFM'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,1,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,0])

#% 3000
ax = fig.add_subplot(436)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_AMJ_fl_PI_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 PI temp AMJ'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,1,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,0])

# 2000
ax = fig.add_subplot(439)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_JAS_fl_PI_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 PI temp JAS'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[1,0,0,0])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,1])    


# 1000
ax = fig.add_subplot(4,3,12)
pos = ax.get_position()
bnd = list(pos.bounds)
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)
cmap = plt.get_cmap('gist_ncar')
ax.set_facecolor('grey')
m.drawcoastlines(linewidth=0.05)
m.fillcontinents(color='white')
contf = m.contourf(x, y, temp_OND_fl_PI_0[0,:,:], contf_lvls, cmap=cmap, extend='both')
title_name = 'KDS75 PI temp OND'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,1])

cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe, orientation='horizontal', drawedges=True)
cbar.set_label('degree C') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)

output_ls = os.listdir(figures_path)

scriptname = 'p02_plot_MMM_TS_maps.py'

#
if scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:2] \
            + '_fig1_' + '.png', bbox_inches='tight', dpi=300)

#plt.close('all') # close







