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
def round_to_base(x, base=5):
    return int(base * round(float(x) / base))
#
import sys
#
import pickle

os.chdir('/Users/earl/SAMexp')
figures_path = '/Users/earl/Dropbox/SAMexp/figures/'
data_path = '/Users/earl/Dropbox/Data/SAMexp/'

MMM = ['JFM', 'AMJ', 'JAS', 'OND']

DATA = ['', '_UP', '_PI']
DATA_out = ['CT', 'UP', 'PI']

scriptname = os.path.basename(sys.argv[0])[:-3]

out = {'out':'out'}

print('OK, all set')


#%% Load data temp_x
z_idx = 0

for d, d_out in zip(DATA, DATA_out):
    for m in MMM:
        # create an instance of the ncCDF4 class
        nc_fid = nc.Dataset(data_path + 'KDS75' + d + '/salt_'
                            + m + '_0.nc', 'r')
        # same for temperature (1,2,3)
        out[d_out + m] = nc_fid.variables['salt'][0,z_idx,:,:]
        
        #
        print(d_out + m + ' OK !')
            
# get dimensions
lat = nc_fid.variables['yt_ocean'][:]
lon1 = nc_fid.variables['xt_ocean'][:]
z = nc_fid.variables['st_ocean'][z_idx]


#%% re-organise longitudes
out1 = {'1':'re-arrange data so that maps start at 70W'}
lon_less_m70 = lon1 < -70
lon_less_m70_flipped = np.flipud(lon_less_m70)

for d, d_out in zip(DATA, DATA_out):
    for m in MMM:
        ma = np.ma.empty((988,3600))
        ma[:,lon_less_m70_flipped] = out[d_out + m][:,lon_less_m70]
        ma[:,~lon_less_m70_flipped] = out[d_out + m][:,~lon_less_m70]
        out1[d_out + m] = ma
        print(d_out + m + ' OK !')


#%% prepare data for plot
outf = {'f':'create anomalies for columns 1 and 2'}

for d, d_out in zip(DATA, DATA_out):
    for m in MMM:
        if d_out is 'CT':
            outf[d_out + m] = out1[d_out + m]
            
        else:
            outf[d_out + m] = out1[d_out + m] - out1['CT' + m]
        
        print(d_out + m + ' OK !')
            
            
#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
Bm = Basemap(projection='mill', llcrnrlat=-55.01,urcrnrlat=-19.99,\
llcrnrlon=99.99,urcrnrlon=170.01, resolution='c')

pickle.dump(Bm,open('map.pickle','wb'),-1)  # pickle it 

Bm.fix_aspect = True
m_aspect = Bm.aspect

lon = lon1 + 210


#%% plot surface maps: TEMPERATURE
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 14}) 

fig.set_size_inches(12, 10) # set figure size in inches
row = 4
col = 3

merid = np.zeros(12)
paral = np.zeros(12)

merid[-col:] = 1
paral[np.arange(0,11,col)] = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

s = 0

seasons =['Summer', 'Autumn', 'Winter', 'Spring']

exps = ['CT', 'UP', 'PI']
ex = 0

for m in MMM:
    for d, d_out in zip(DATA, DATA_out):
        s = s + 1

        # position figure wrt window template
        ax = fig.add_subplot(row,col,s)
        pos = ax.get_position()
        bnd = list(pos.bounds)
        magn = 0.04
        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
        ax.set_position(bnd)
        
        # colourmap: blue to red because looking at bias.
        if d_out is 'CT':
            cmap = plt.get_cmap('gist_ncar')
            step = 0.2
            # levels to show on colourbar
            contf_lvls = np.arange(33.4,37+1e-08,step)
            
        else:
            cmap = plt.get_cmap('seismic')
            step = 0.2
            # levels to show on colourbar
            contf_lvls = np.arange(-1.6,1.6+1e-08,step)
        
        
        ax.set_facecolor('grey')
        
        # meshgrid of lon lats
        lons, lats = np.meshgrid(lon, lat)
        # use projection template to create x and y axis
        x, y = Bm(lons, lats)
        
        #
        pickle.load(open('map.pickle','rb'))   # load here the above pickle
        
        # draw land outlines
        Bm.drawcoastlines(linewidth=0.05)
        Bm.fillcontinents(color='white')
        
        # filled contour plot
        contf = Bm.contourf(x, y, outf[d_out + m], 
                           contf_lvls, cmap=cmap, extend='both')
        
        # title ...
        if s <= 3:
            title_name = exps[s-1]
            ax.set_title(title_name)
            
        if s is 1 or s is 4 or s is 7 or s is 10:
            ex = ex + 1
            ax.set_ylabel(seasons[ex-1])
            ax.yaxis.labelpad = 35
        
        # meridians. last input is meridians tick label
        meridians = map(round_to_base, np.arange(110, 180, 20))
        Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,merid[s-1]])
        # parallels. last input is paralles tick label
        parallels = map(round_to_base, np.arange(-50, -10, 10))
        Bm.drawparallels(parallels, linewidth=0.2, labels=[paral[s-1],0,0,0])
        
        if s is 10:
            cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*4] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe,
                                orientation='vertical', drawedges=True)
            cbar.set_label(r'Absolute S $(psu)$') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.yaxis.set_ticks_position('left')
            
        elif s is 12:
            cbar_pos = [bnd[0]+bnd[2]+0.01, bnd[1], 0.01, bnd[3]*4] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='vertical', drawedges=True)
            cbar.set_label(r'S anomaly $(psu)$') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        print(d_out + m + ' OK !')
        
plt.suptitle(r"Sea Surface Salinity", y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'
    
#
if scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig1_' + '.png', bbox_inches='tight', dpi=100)

