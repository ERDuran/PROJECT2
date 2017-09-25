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


#% Load data tau_x
# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Jan' + '.nc', 'r')

# same for temperature (1,2,3)
tau_x_Jan = nc_fid.variables['tau_x'][0,:,:]

# Feb
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Feb' + '.nc', 'r')
tau_x_Feb = nc_fid.variables['tau_x'][0,:,:]
# Mar
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Mar' + '.nc', 'r')
tau_x_Mar = nc_fid.variables['tau_x'][0,:,:]

# MMM average
tau_x_JFM_CT = (tau_x_Jan + tau_x_Feb + tau_x_Mar)/3


# AMJ
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Apr' + '.nc', 'r')
tau_x_Apr = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_May' + '.nc', 'r')
tau_x_May = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Jun' + '.nc', 'r')
tau_x_Jun = nc_fid.variables['tau_x'][0,:,:]
tau_x_AMJ_CT = (tau_x_Apr + tau_x_May + tau_x_Jun)/3


# JAS
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Jul' + '.nc', 'r')
tau_x_Jul = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Aug' + '.nc', 'r')
tau_x_Aug = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Sep' + '.nc', 'r')
tau_x_Sep = nc_fid.variables['tau_x'][0,:,:]
tau_x_JAS_CT = (tau_x_Jul + tau_x_Aug + tau_x_Sep)/3


# OND
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Oct' + '.nc', 'r')
tau_x_Oct = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Nov' + '.nc', 'r')
tau_x_Nov = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_x_Dec' + '.nc', 'r')
tau_x_Dec = nc_fid.variables['tau_x'][0,:,:]
tau_x_OND_CT = (tau_x_Oct + tau_x_Nov + tau_x_Dec)/3

# get dimensions
lat = nc_fid.variables['yu_ocean'][:]
lon1 = nc_fid.variables['xu_ocean'][:]


#% Load data tau_y
# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Jan' + '.nc', 'r')

# same for temperature (1,2,3)
tau_y_Jan = nc_fid.variables['tau_y'][0,:,:]

# Feb
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Feb' + '.nc', 'r')
tau_y_Feb = nc_fid.variables['tau_y'][0,:,:]
# Mar
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Mar' + '.nc', 'r')
tau_y_Mar = nc_fid.variables['tau_y'][0,:,:]

# MMM average
tau_y_JFM_CT = (tau_y_Jan + tau_y_Feb + tau_y_Mar)/3


# AMJ
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Apr' + '.nc', 'r')
tau_y_Apr = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_May' + '.nc', 'r')
tau_y_May = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Jun' + '.nc', 'r')
tau_y_Jun = nc_fid.variables['tau_y'][0,:,:]
tau_y_AMJ_CT = (tau_y_Apr + tau_y_May + tau_y_Jun)/3


# JAS
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Jul' + '.nc', 'r')
tau_y_Jul = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Aug' + '.nc', 'r')
tau_y_Aug = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Sep' + '.nc', 'r')
tau_y_Sep = nc_fid.variables['tau_y'][0,:,:]
tau_y_JAS_CT = (tau_y_Jul + tau_y_Aug + tau_y_Sep)/3


# OND
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Oct' + '.nc', 'r')
tau_y_Oct = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Nov' + '.nc', 'r')
tau_y_Nov = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75/tau_y_Dec' + '.nc', 'r')
tau_y_Dec = nc_fid.variables['tau_y'][0,:,:]
tau_y_OND_CT = (tau_y_Oct + tau_y_Nov + tau_y_Dec)/3


#% UP
# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Jan' + '.nc', 'r')
tau_x_Jan = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Feb' + '.nc', 'r')
tau_x_Feb = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Mar' + '.nc', 'r')
tau_x_Mar = nc_fid.variables['tau_x'][0,:,:]
tau_x_JFM_UP = (tau_x_Jan + tau_x_Feb + tau_x_Mar)/3


# AMJ
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Apr' + '.nc', 'r')
tau_x_Apr = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_May' + '.nc', 'r')
tau_x_May = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Jun' + '.nc', 'r')
tau_x_Jun = nc_fid.variables['tau_x'][0,:,:]
tau_x_AMJ_UP = (tau_x_Apr + tau_x_May + tau_x_Jun)/3


# JAS
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Jul' + '.nc', 'r')
tau_x_Jul = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Aug' + '.nc', 'r')
tau_x_Aug = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Sep' + '.nc', 'r')
tau_x_Sep = nc_fid.variables['tau_x'][0,:,:]
tau_x_JAS_UP = (tau_x_Jul + tau_x_Aug + tau_x_Sep)/3


# OND
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Oct' + '.nc', 'r')
tau_x_Oct = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Nov' + '.nc', 'r')
tau_x_Nov = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_x_Dec' + '.nc', 'r')
tau_x_Dec = nc_fid.variables['tau_x'][0,:,:]
tau_x_OND_UP = (tau_x_Oct + tau_x_Nov + tau_x_Dec)/3


# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Jan' + '.nc', 'r')
tau_y_Jan = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Feb' + '.nc', 'r')
tau_y_Feb = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Mar' + '.nc', 'r')
tau_y_Mar = nc_fid.variables['tau_y'][0,:,:]
tau_y_JFM_UP = (tau_y_Jan + tau_y_Feb + tau_y_Mar)/3


# AMJ
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Apr' + '.nc', 'r')
tau_y_Apr = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_May' + '.nc', 'r')
tau_y_May = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Jun' + '.nc', 'r')
tau_y_Jun = nc_fid.variables['tau_y'][0,:,:]
tau_y_AMJ_UP = (tau_y_Apr + tau_y_May + tau_y_Jun)/3


# JAS
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Jul' + '.nc', 'r')
tau_y_Jul = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Aug' + '.nc', 'r')
tau_y_Aug = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Sep' + '.nc', 'r')
tau_y_Sep = nc_fid.variables['tau_y'][0,:,:]
tau_y_JAS_UP = (tau_y_Jul + tau_y_Aug + tau_y_Sep)/3


# OND
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Oct' + '.nc', 'r')
tau_y_Oct = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Nov' + '.nc', 'r')
tau_y_Nov = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_UP/tau_y_Dec' + '.nc', 'r')
tau_y_Dec = nc_fid.variables['tau_y'][0,:,:]
tau_y_OND_UP = (tau_y_Oct + tau_y_Nov + tau_y_Dec)/3


#% PI
# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Jan' + '.nc', 'r')
tau_x_Jan = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Feb' + '.nc', 'r')
tau_x_Feb = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Mar' + '.nc', 'r')
tau_x_Mar = nc_fid.variables['tau_x'][0,:,:]
tau_x_JFM_PI = (tau_x_Jan + tau_x_Feb + tau_x_Mar)/3


# AMJ
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Apr' + '.nc', 'r')
tau_x_Apr = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_May' + '.nc', 'r')
tau_x_May = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Jun' + '.nc', 'r')
tau_x_Jun = nc_fid.variables['tau_x'][0,:,:]
tau_x_AMJ_PI = (tau_x_Apr + tau_x_May + tau_x_Jun)/3


# JAS
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Jul' + '.nc', 'r')
tau_x_Jul = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Aug' + '.nc', 'r')
tau_x_Aug = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Sep' + '.nc', 'r')
tau_x_Sep = nc_fid.variables['tau_x'][0,:,:]
tau_x_JAS_PI = (tau_x_Jul + tau_x_Aug + tau_x_Sep)/3


# OND
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Oct' + '.nc', 'r')
tau_x_Oct = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Nov' + '.nc', 'r')
tau_x_Nov = nc_fid.variables['tau_x'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_x_Dec' + '.nc', 'r')
tau_x_Dec = nc_fid.variables['tau_x'][0,:,:]
tau_x_OND_PI = (tau_x_Oct + tau_x_Nov + tau_x_Dec)/3


# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Jan' + '.nc', 'r')
tau_y_Jan = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Feb' + '.nc', 'r')
tau_y_Feb = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Mar' + '.nc', 'r')
tau_y_Mar = nc_fid.variables['tau_y'][0,:,:]
tau_y_JFM_PI = (tau_y_Jan + tau_y_Feb + tau_y_Mar)/3


# AMJ
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Apr' + '.nc', 'r')
tau_y_Apr = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_May' + '.nc', 'r')
tau_y_May = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Jun' + '.nc', 'r')
tau_y_Jun = nc_fid.variables['tau_y'][0,:,:]
tau_y_AMJ_PI = (tau_y_Apr + tau_y_May + tau_y_Jun)/3


# JAS
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Jul' + '.nc', 'r')
tau_y_Jul = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Aug' + '.nc', 'r')
tau_y_Aug = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Sep' + '.nc', 'r')
tau_y_Sep = nc_fid.variables['tau_y'][0,:,:]
tau_y_JAS_PI = (tau_y_Jul + tau_y_Aug + tau_y_Sep)/3


# OND
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Oct' + '.nc', 'r')
tau_y_Oct = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Nov' + '.nc', 'r')
tau_y_Nov = nc_fid.variables['tau_y'][0,:,:]
nc_fid = nc.Dataset(data_path + 'KDS75_PI/tau_y_Dec' + '.nc', 'r')
tau_y_Dec = nc_fid.variables['tau_y'][0,:,:]
tau_y_OND_PI = (tau_y_Oct + tau_y_Nov + tau_y_Dec)/3



#%% re-organise longitudes
lon_less_m180 = lon1 < -180

lon_less_m180_flipped = np.flipud(lon_less_m180)

tau_x_JFM_fl1_CT = np.ma.empty((989,3600))
tau_x_JFM_fl1_CT[:,lon_less_m180_flipped] = tau_x_JFM_CT[:,lon_less_m180];
tau_x_JFM_fl1_CT[:,~lon_less_m180_flipped] = tau_x_JFM_CT[:,~lon_less_m180];

tau_x_AMJ_fl1_CT = np.ma.empty((989,3600))
tau_x_AMJ_fl1_CT[:,lon_less_m180_flipped] = tau_x_AMJ_CT[:,lon_less_m180];
tau_x_AMJ_fl1_CT[:,~lon_less_m180_flipped] = tau_x_AMJ_CT[:,~lon_less_m180];

tau_x_JAS_fl1_CT = np.ma.empty((989,3600))
tau_x_JAS_fl1_CT[:,lon_less_m180_flipped] = tau_x_JAS_CT[:,lon_less_m180];
tau_x_JAS_fl1_CT[:,~lon_less_m180_flipped] = tau_x_JAS_CT[:,~lon_less_m180];

tau_x_OND_fl1_CT = np.ma.empty((989,3600))
tau_x_OND_fl1_CT[:,lon_less_m180_flipped] = tau_x_OND_CT[:,lon_less_m180];
tau_x_OND_fl1_CT[:,~lon_less_m180_flipped] = tau_x_OND_CT[:,~lon_less_m180];


#%
tau_x_JFM_fl1_UP = np.ma.empty((989,3600))
tau_x_JFM_fl1_UP[:,lon_less_m180_flipped] = tau_x_JFM_UP[:,lon_less_m180];
tau_x_JFM_fl1_UP[:,~lon_less_m180_flipped] = tau_x_JFM_UP[:,~lon_less_m180];

tau_x_AMJ_fl1_UP = np.ma.empty((989,3600))
tau_x_AMJ_fl1_UP[:,lon_less_m180_flipped] = tau_x_AMJ_UP[:,lon_less_m180];
tau_x_AMJ_fl1_UP[:,~lon_less_m180_flipped] = tau_x_AMJ_UP[:,~lon_less_m180];

tau_x_JAS_fl1_UP = np.ma.empty((989,3600))
tau_x_JAS_fl1_UP[:,lon_less_m180_flipped] = tau_x_JAS_UP[:,lon_less_m180];
tau_x_JAS_fl1_UP[:,~lon_less_m180_flipped] = tau_x_JAS_UP[:,~lon_less_m180];

tau_x_OND_fl1_UP = np.ma.empty((989,3600))
tau_x_OND_fl1_UP[:,lon_less_m180_flipped] = tau_x_OND_UP[:,lon_less_m180];
tau_x_OND_fl1_UP[:,~lon_less_m180_flipped] = tau_x_OND_UP[:,~lon_less_m180];


#%
tau_x_JFM_fl1_PI = np.ma.empty((989,3600))
tau_x_JFM_fl1_PI[:,lon_less_m180_flipped] = tau_x_JFM_PI[:,lon_less_m180];
tau_x_JFM_fl1_PI[:,~lon_less_m180_flipped] = tau_x_JFM_PI[:,~lon_less_m180];

tau_x_AMJ_fl1_PI = np.ma.empty((989,3600))
tau_x_AMJ_fl1_PI[:,lon_less_m180_flipped] = tau_x_AMJ_PI[:,lon_less_m180];
tau_x_AMJ_fl1_PI[:,~lon_less_m180_flipped] = tau_x_AMJ_PI[:,~lon_less_m180];

tau_x_JAS_fl1_PI = np.ma.empty((989,3600))
tau_x_JAS_fl1_PI[:,lon_less_m180_flipped] = tau_x_JAS_PI[:,lon_less_m180];
tau_x_JAS_fl1_PI[:,~lon_less_m180_flipped] = tau_x_JAS_PI[:,~lon_less_m180];

tau_x_OND_fl1_PI = np.ma.empty((989,3600))
tau_x_OND_fl1_PI[:,lon_less_m180_flipped] = tau_x_OND_PI[:,lon_less_m180];
tau_x_OND_fl1_PI[:,~lon_less_m180_flipped] = tau_x_OND_PI[:,~lon_less_m180];



#%% re-organise longitudes
lon2 = lon1 + 100
lon_less_m70 = lon2 < -70

lon_less_m70_flipped = np.flipud(lon_less_m70)

tau_x_JFM_fl_CT = np.ma.empty((989,3600))
tau_x_JFM_fl_CT[:,lon_less_m70_flipped] = tau_x_JFM_fl1_CT[:,lon_less_m70];
tau_x_JFM_fl_CT[:,~lon_less_m70_flipped] = tau_x_JFM_fl1_CT[:,~lon_less_m70];

tau_x_AMJ_fl_CT = np.ma.empty((989,3600))
tau_x_AMJ_fl_CT[:,lon_less_m70_flipped] = tau_x_AMJ_fl1_CT[:,lon_less_m70];
tau_x_AMJ_fl_CT[:,~lon_less_m70_flipped] = tau_x_AMJ_fl1_CT[:,~lon_less_m70];

tau_x_JAS_fl_CT = np.ma.empty((989,3600))
tau_x_JAS_fl_CT[:,lon_less_m70_flipped] = tau_x_JAS_fl1_CT[:,lon_less_m70];
tau_x_JAS_fl_CT[:,~lon_less_m70_flipped] = tau_x_JAS_fl1_CT[:,~lon_less_m70];

tau_x_OND_fl_CT = np.ma.empty((989,3600))
tau_x_OND_fl_CT[:,lon_less_m70_flipped] = tau_x_OND_fl1_CT[:,lon_less_m70];
tau_x_OND_fl_CT[:,~lon_less_m70_flipped] = tau_x_OND_fl1_CT[:,~lon_less_m70];


#%
tau_x_JFM_fl_UP = np.ma.empty((989,3600))
tau_x_JFM_fl_UP[:,lon_less_m70_flipped] = tau_x_JFM_fl1_UP[:,lon_less_m70];
tau_x_JFM_fl_UP[:,~lon_less_m70_flipped] = tau_x_JFM_fl1_UP[:,~lon_less_m70];

tau_x_AMJ_fl_UP = np.ma.empty((989,3600))
tau_x_AMJ_fl_UP[:,lon_less_m70_flipped] = tau_x_AMJ_fl1_UP[:,lon_less_m70];
tau_x_AMJ_fl_UP[:,~lon_less_m70_flipped] = tau_x_AMJ_fl1_UP[:,~lon_less_m70];

tau_x_JAS_fl_UP = np.ma.empty((989,3600))
tau_x_JAS_fl_UP[:,lon_less_m70_flipped] = tau_x_JAS_fl1_UP[:,lon_less_m70];
tau_x_JAS_fl_UP[:,~lon_less_m70_flipped] = tau_x_JAS_fl1_UP[:,~lon_less_m70];

tau_x_OND_fl_UP = np.ma.empty((989,3600))
tau_x_OND_fl_UP[:,lon_less_m70_flipped] = tau_x_OND_fl1_UP[:,lon_less_m70];
tau_x_OND_fl_UP[:,~lon_less_m70_flipped] = tau_x_OND_fl1_UP[:,~lon_less_m70];


#%
tau_x_JFM_fl_PI = np.ma.empty((989,3600))
tau_x_JFM_fl_PI[:,lon_less_m70_flipped] = tau_x_JFM_fl1_PI[:,lon_less_m70];
tau_x_JFM_fl_PI[:,~lon_less_m70_flipped] = tau_x_JFM_fl1_PI[:,~lon_less_m70];

tau_x_AMJ_fl_PI = np.ma.empty((989,3600))
tau_x_AMJ_fl_PI[:,lon_less_m70_flipped] = tau_x_AMJ_fl1_PI[:,lon_less_m70];
tau_x_AMJ_fl_PI[:,~lon_less_m70_flipped] = tau_x_AMJ_fl1_PI[:,~lon_less_m70];

tau_x_JAS_fl_PI = np.ma.empty((989,3600))
tau_x_JAS_fl_PI[:,lon_less_m70_flipped] = tau_x_JAS_fl1_PI[:,lon_less_m70];
tau_x_JAS_fl_PI[:,~lon_less_m70_flipped] = tau_x_JAS_fl1_PI[:,~lon_less_m70];

tau_x_OND_fl_PI = np.ma.empty((989,3600))
tau_x_OND_fl_PI[:,lon_less_m70_flipped] = tau_x_OND_fl1_PI[:,lon_less_m70];
tau_x_OND_fl_PI[:,~lon_less_m70_flipped] = tau_x_OND_fl1_PI[:,~lon_less_m70];



#%% re-organise TAU_Y
tau_y_JFM_fl1_CT = np.ma.empty((989,3600))
tau_y_JFM_fl1_CT[:,lon_less_m180_flipped] = tau_y_JFM_CT[:,lon_less_m180];
tau_y_JFM_fl1_CT[:,~lon_less_m180_flipped] = tau_y_JFM_CT[:,~lon_less_m180];
tau_y_AMJ_fl1_CT = np.ma.empty((989,3600))
tau_y_AMJ_fl1_CT[:,lon_less_m180_flipped] = tau_y_AMJ_CT[:,lon_less_m180];
tau_y_AMJ_fl1_CT[:,~lon_less_m180_flipped] = tau_y_AMJ_CT[:,~lon_less_m180];
tau_y_JAS_fl1_CT = np.ma.empty((989,3600))
tau_y_JAS_fl1_CT[:,lon_less_m180_flipped] = tau_y_JAS_CT[:,lon_less_m180];
tau_y_JAS_fl1_CT[:,~lon_less_m180_flipped] = tau_y_JAS_CT[:,~lon_less_m180];
tau_y_OND_fl1_CT = np.ma.empty((989,3600))
tau_y_OND_fl1_CT[:,lon_less_m180_flipped] = tau_y_OND_CT[:,lon_less_m180];
tau_y_OND_fl1_CT[:,~lon_less_m180_flipped] = tau_y_OND_CT[:,~lon_less_m180];
tau_y_JFM_fl1_UP = np.ma.empty((989,3600))
tau_y_JFM_fl1_UP[:,lon_less_m180_flipped] = tau_y_JFM_UP[:,lon_less_m180];
tau_y_JFM_fl1_UP[:,~lon_less_m180_flipped] = tau_y_JFM_UP[:,~lon_less_m180];
tau_y_AMJ_fl1_UP = np.ma.empty((989,3600))
tau_y_AMJ_fl1_UP[:,lon_less_m180_flipped] = tau_y_AMJ_UP[:,lon_less_m180];
tau_y_AMJ_fl1_UP[:,~lon_less_m180_flipped] = tau_y_AMJ_UP[:,~lon_less_m180];
tau_y_JAS_fl1_UP = np.ma.empty((989,3600))
tau_y_JAS_fl1_UP[:,lon_less_m180_flipped] = tau_y_JAS_UP[:,lon_less_m180];
tau_y_JAS_fl1_UP[:,~lon_less_m180_flipped] = tau_y_JAS_UP[:,~lon_less_m180];
tau_y_OND_fl1_UP = np.ma.empty((989,3600))
tau_y_OND_fl1_UP[:,lon_less_m180_flipped] = tau_y_OND_UP[:,lon_less_m180];
tau_y_OND_fl1_UP[:,~lon_less_m180_flipped] = tau_y_OND_UP[:,~lon_less_m180];
tau_y_JFM_fl1_PI = np.ma.empty((989,3600))
tau_y_JFM_fl1_PI[:,lon_less_m180_flipped] = tau_y_JFM_PI[:,lon_less_m180];
tau_y_JFM_fl1_PI[:,~lon_less_m180_flipped] = tau_y_JFM_PI[:,~lon_less_m180];
tau_y_AMJ_fl1_PI = np.ma.empty((989,3600))
tau_y_AMJ_fl1_PI[:,lon_less_m180_flipped] = tau_y_AMJ_PI[:,lon_less_m180];
tau_y_AMJ_fl1_PI[:,~lon_less_m180_flipped] = tau_y_AMJ_PI[:,~lon_less_m180];
tau_y_JAS_fl1_PI = np.ma.empty((989,3600))
tau_y_JAS_fl1_PI[:,lon_less_m180_flipped] = tau_y_JAS_PI[:,lon_less_m180];
tau_y_JAS_fl1_PI[:,~lon_less_m180_flipped] = tau_y_JAS_PI[:,~lon_less_m180];
tau_y_OND_fl1_PI = np.ma.empty((989,3600))
tau_y_OND_fl1_PI[:,lon_less_m180_flipped] = tau_y_OND_PI[:,lon_less_m180];
tau_y_OND_fl1_PI[:,~lon_less_m180_flipped] = tau_y_OND_PI[:,~lon_less_m180];
tau_y_JFM_fl_CT = np.ma.empty((989,3600))
tau_y_JFM_fl_CT[:,lon_less_m70_flipped] = tau_y_JFM_fl1_CT[:,lon_less_m70];
tau_y_JFM_fl_CT[:,~lon_less_m70_flipped] = tau_y_JFM_fl1_CT[:,~lon_less_m70];
tau_y_AMJ_fl_CT = np.ma.empty((989,3600))
tau_y_AMJ_fl_CT[:,lon_less_m70_flipped] = tau_y_AMJ_fl1_CT[:,lon_less_m70];
tau_y_AMJ_fl_CT[:,~lon_less_m70_flipped] = tau_y_AMJ_fl1_CT[:,~lon_less_m70];
tau_y_JAS_fl_CT = np.ma.empty((989,3600))
tau_y_JAS_fl_CT[:,lon_less_m70_flipped] = tau_y_JAS_fl1_CT[:,lon_less_m70];
tau_y_JAS_fl_CT[:,~lon_less_m70_flipped] = tau_y_JAS_fl1_CT[:,~lon_less_m70];
tau_y_OND_fl_CT = np.ma.empty((989,3600))
tau_y_OND_fl_CT[:,lon_less_m70_flipped] = tau_y_OND_fl1_CT[:,lon_less_m70];
tau_y_OND_fl_CT[:,~lon_less_m70_flipped] = tau_y_OND_fl1_CT[:,~lon_less_m70];
tau_y_JFM_fl_UP = np.ma.empty((989,3600))
tau_y_JFM_fl_UP[:,lon_less_m70_flipped] = tau_y_JFM_fl1_UP[:,lon_less_m70];
tau_y_JFM_fl_UP[:,~lon_less_m70_flipped] = tau_y_JFM_fl1_UP[:,~lon_less_m70];
tau_y_AMJ_fl_UP = np.ma.empty((989,3600))
tau_y_AMJ_fl_UP[:,lon_less_m70_flipped] = tau_y_AMJ_fl1_UP[:,lon_less_m70];
tau_y_AMJ_fl_UP[:,~lon_less_m70_flipped] = tau_y_AMJ_fl1_UP[:,~lon_less_m70];
tau_y_JAS_fl_UP = np.ma.empty((989,3600))
tau_y_JAS_fl_UP[:,lon_less_m70_flipped] = tau_y_JAS_fl1_UP[:,lon_less_m70];
tau_y_JAS_fl_UP[:,~lon_less_m70_flipped] = tau_y_JAS_fl1_UP[:,~lon_less_m70];
tau_y_OND_fl_UP = np.ma.empty((989,3600))
tau_y_OND_fl_UP[:,lon_less_m70_flipped] = tau_y_OND_fl1_UP[:,lon_less_m70];
tau_y_OND_fl_UP[:,~lon_less_m70_flipped] = tau_y_OND_fl1_UP[:,~lon_less_m70];
tau_y_JFM_fl_PI = np.ma.empty((989,3600))
tau_y_JFM_fl_PI[:,lon_less_m70_flipped] = tau_y_JFM_fl1_PI[:,lon_less_m70];
tau_y_JFM_fl_PI[:,~lon_less_m70_flipped] = tau_y_JFM_fl1_PI[:,~lon_less_m70];
tau_y_AMJ_fl_PI = np.ma.empty((989,3600))
tau_y_AMJ_fl_PI[:,lon_less_m70_flipped] = tau_y_AMJ_fl1_PI[:,lon_less_m70];
tau_y_AMJ_fl_PI[:,~lon_less_m70_flipped] = tau_y_AMJ_fl1_PI[:,~lon_less_m70];
tau_y_JAS_fl_PI = np.ma.empty((989,3600))
tau_y_JAS_fl_PI[:,lon_less_m70_flipped] = tau_y_JAS_fl1_PI[:,lon_less_m70];
tau_y_JAS_fl_PI[:,~lon_less_m70_flipped] = tau_y_JAS_fl1_PI[:,~lon_less_m70];
tau_y_OND_fl_PI = np.ma.empty((989,3600))
tau_y_OND_fl_PI[:,lon_less_m70_flipped] = tau_y_OND_fl1_PI[:,lon_less_m70];
tau_y_OND_fl_PI[:,~lon_less_m70_flipped] = tau_y_OND_fl1_PI[:,~lon_less_m70];



#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
m = Basemap(projection='gall', llcrnrlat=-70,urcrnrlat=0.01,\
llcrnrlon=-70,urcrnrlon=290, resolution='h')
m.fix_aspect = True
m_aspect = m.aspect

lon = lon2 + 110


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
contf_lvls = np.arange(-0.10,0.405,0.05) # levels to show on colourbar
contf = m.contourf(x, y, tau_x_JFM_fl_CT, contf_lvls, cmap=cmap, extend='both')

# quiver
m.quiver(x[points], y[points], 
    tau_x_JFM_fl_CT[points], tau_y_JFM_fl_CT[points], scale=quiv_scale)

# title ...
title_name = 'KDS75 tau x JFM'
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
contf = m.contourf(x, y, tau_x_AMJ_fl_CT, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_AMJ_fl_CT[points], tau_y_AMJ_fl_CT[points], scale=quiv_scale)
title_name = 'KDS75 tau x AMJ'
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
contf = m.contourf(x, y, tau_x_JAS_fl_CT, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_JAS_fl_CT[points], tau_y_JAS_fl_CT[points], scale=quiv_scale)
title_name = 'KDS75 tau x JAS'
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
contf = m.contourf(x, y, tau_x_OND_fl_CT, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_OND_fl_CT[points], tau_y_OND_fl_CT[points], scale=quiv_scale)
title_name = 'KDS75 tau x OND'
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
contf = m.contourf(x, y, tau_x_JFM_fl_UP, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_JFM_fl_UP[points], tau_y_JFM_fl_UP[points], scale=quiv_scale)
title_name = 'KDS75 UP tau x JFM'
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
contf = m.contourf(x, y, tau_x_AMJ_fl_UP, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_AMJ_fl_UP[points], tau_y_AMJ_fl_UP[points], scale=quiv_scale)
title_name = 'KDS75 UP tau x AMJ'
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
contf = m.contourf(x, y, tau_x_JAS_fl_UP, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_JAS_fl_UP[points], tau_y_JAS_fl_UP[points], scale=quiv_scale)
title_name = 'KDS75 UP tau x JAS'
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
contf = m.contourf(x, y, tau_x_OND_fl_UP, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_OND_fl_UP[points], tau_y_OND_fl_UP[points], scale=quiv_scale)
title_name = 'KDS75 UP tau x OND'
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
contf = m.contourf(x, y, tau_x_JFM_fl_PI, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_JFM_fl_PI[points], tau_y_JFM_fl_PI[points], scale=quiv_scale)
title_name = 'KDS75 PI tau x JFM'
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
contf = m.contourf(x, y, tau_x_AMJ_fl_PI, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_AMJ_fl_PI[points], tau_y_AMJ_fl_PI[points], scale=quiv_scale)
title_name = 'KDS75 PI tau x AMJ'
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
contf = m.contourf(x, y, tau_x_JAS_fl_PI, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_JAS_fl_PI[points], tau_y_JAS_fl_PI[points], scale=quiv_scale)
title_name = 'KDS75 PI tau x JAS'
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
contf = m.contourf(x, y, tau_x_OND_fl_PI, contf_lvls, cmap=cmap, extend='both')
m.quiver(x[points], y[points], 
    tau_x_OND_fl_PI[points], tau_y_OND_fl_PI[points], scale=quiv_scale)
title_name = 'KDS75 PI tau x OND'
ax.set_title(title_name, y=0.97)
meridians = map(round_to_base, np.arange(-160, 180, 40))
m.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
parallels = map(round_to_base, np.arange(-80, 100, 20))
m.drawparallels(parallels, linewidth=0.2, labels=[0,0,0,1])

cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe, orientation='horizontal', drawedges=True)
cbar.set_label('N/m^2') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)

output_ls = os.listdir(figures_path)

#
if scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:2] \
            + '_fig1_' + '.png', bbox_inches='tight', dpi=300)

#plt.close('all') # close







