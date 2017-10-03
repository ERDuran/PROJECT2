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
#
import pickle

os.chdir('/Users/earl/SAMexp')
figures_path = '/Users/earl/Dropbox/SAMexp/figures/'
data_path = '/Users/earl/Dropbox/Data/SAMexp/'

Mon = {'JFM':['Jan', 'Feb', 'Mar'], 'AMJ':['Apr', 'May', 'Jun'], 
       'JAS':['Jul', 'Aug', 'Sep'], 'OND':['Oct', 'Nov', 'Dec']}
 
MMM = ['JFM', 'AMJ', 'JAS', 'OND']

DATA = ['', '_UP', '_PI']
DATA_out = ['CT', 'UP', 'PI']

var = ['x', 'y']

scriptname = os.path.basename(sys.argv[0])[:-3]

temp0 = {'0':'temp at z ranges'}

temp1 = {'1':'temp along latitude sections'}

print('OK, all set')


#%%
# create an instance of the ncCDF4 class
nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_0' + '.nc', 'r')

# get dimensions
lat = nc_fid.variables['yt_ocean'][:]
lon1 = nc_fid.variables['xt_ocean'][:]
lon = lon1 + 210

lon_lines_idx = np.arange(300,3600,400)
lon_lines = lon[lon_lines_idx]

lon_lines_m210 = lon_lines - 210

z0 = nc_fid.variables['st_ocean'][:]

nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_1' + '.nc', 'r')
z1 = nc_fid.variables['st_ocean'][:]

nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_2' + '.nc', 'r')
z2 = nc_fid.variables['st_ocean'][:]

nc_fid = nc.Dataset(data_path + 'KDS75/temp_JFM_3' + '.nc', 'r')
z3 = nc_fid.variables['st_ocean'][:]

z = np.concatenate([z0, z1, z2, z3])

lon_lines = np.round(lon[lon_lines_idx])


#%% Load data tau_x CT
z_idx = 4

for d, d_out in zip(DATA, DATA_out):
    for m in MMM:
        for l,l_out in zip(lon_lines_idx,lon_lines):
            for s in range(0,4):
                # create an instance of the ncCDF4 class
                nc_fid = nc.Dataset(data_path + 'KDS75' + d + 
                                    '/temp_' + m + '_' + str(s) + '.nc', 'r')
                
                # same for temperature (1,2,3)
                temp0['z' + str(s)] = nc_fid.variables['temp'][0,:,:,l]
                
            #
            temp1[d_out + m + str(l_out)] = np.ma.concatenate([\
                  temp0['z0'], temp0['z1'], temp0['z2'], temp0['z3']])
            
            #
            print(d_out + m + str(l_out) + ' OK !')


# get dimensions
lat = nc_fid.variables['yt_ocean'][:]
lon1 = nc_fid.variables['xt_ocean'][lon_lines_idx]
lon = lon1 + 210
            

#%% prepare data for plot
tempf = {'f':'create anomalies for columns 1 and 2'}

for d, d_out in zip(DATA, DATA_out):
    for m in MMM:
        for l_out in lon_lines:
            if d_out is 'CT':
                tempf[d_out + m + str(l_out)] = temp1[d_out + m + str(l_out)]
                
            else:
                tempf[d_out + m + str(l_out)] = \
                temp1[d_out + m + str(l_out)] - \
                temp1['CT' + m + str(l_out)]
            
            print(d_out + m + str(l_out) + ' OK !')
            

#%% plot surface maps: TEMPERATURE
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 6}) 

fig.set_size_inches(10.25, 5) # set figure size in inches
row = 1
col = 1

s = 0

for m in MMM[0:1]:
    for d, d_out in zip(DATA[0:1], DATA_out[0:1]):
        for l_out in lon_lines[0:1]:
            
            s = s + 1
            
            # position figure wrt window template
            ax = fig.add_subplot(row,col,s)
    #        pos = ax.get_position()
    #        bnd = list(pos.bounds)
    #        magn = 0.04
    ##        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn]
    #        ax.set_position(bnd)
            
            # colourmap: blue to red because looking at bias.
            if d_out is 'CT':
                cmap = plt.get_cmap('gist_ncar')
                step = 2
                # levels to show on colourbar
                contf_lvls = np.arange(-2,30+1e-08,step)
                
            else:
                cmap = plt.get_cmap('seismic')
                step = 0.3
                # levels to show on colourbar
                contf_lvls = np.arange(-3,3+1e-08,step)            
            
            
            ax.set_facecolor('grey')
            
            # meshgrid of lon lats
            lons, lats = np.meshgrid(lon, lat)
    
            # filled contour plot
            contf = plt.contourf(lat, z, tempf[d_out + m + str(l_out)], 
                               contf_lvls, cmap=cmap, extend='both')
            
            # title ...
            title_name = 'KDS75 ' + d_out + ' SST ' + m + ' z = ' + str(z)
            ax.set_title(title_name)
            
    #        if s is 10:
    ##            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
    #            cbar_axe = fig.add_axes(cbar_pos)
    #            cbar = plt.colorbar(contf, cax=cbar_axe, 
    #                                orientation='horizontal', drawedges=True)
    #            cbar.set_label('degree C') # units label on colourbar
    #            cbar.dividers.set_linewidth(0.2)
    #            cbar.outline.set_linewidth(0.5)
    #            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
    #            cbar.ax.tick_params(width=0.2, length= 2)
                
    #        elif s is 11:
    ##            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
    #            cbar_axe = fig.add_axes(cbar_pos)
    #            cbar = plt.colorbar(contf, cax=cbar_axe, 
    #                                orientation='horizontal', drawedges=True)
    #            cbar.set_label('degree C') # units label on colourbar
    #            cbar.dividers.set_linewidth(0.2)
    #            cbar.outline.set_linewidth(0.5)
    #            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
    #            cbar.ax.tick_params(width=0.2, length= 2)
                
            print(d_out + m + ' OK !')
        

output_ls = os.listdir(figures_path)


#
if scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)

# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig1_' + str(z_idx) + '.png', bbox_inches='tight', dpi=300)

