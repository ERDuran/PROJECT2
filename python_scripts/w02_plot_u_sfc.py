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

var = ['u', 'v']

scriptname = os.path.basename(sys.argv[0])[:-3]

out = {'out':'out'}

print('OK, all set')


#%% Load data temp_x
z_idx = 0

for v in var:
    for d, d_out in zip(DATA, DATA_out):
        for m in MMM:
            # create an instance of the ncCDF4 class
            nc_fid = nc.Dataset(data_path + 'KDS75' + d + '/' + v + '_'
                                + m + '_0.nc', 'r')
            # same for temperature (1,2,3)
            out[v + d_out + m] = nc_fid.variables[v][0,z_idx,:,:]
            
            #
            print(v + d_out + m + ' OK !')
            
# get dimensions
lat = nc_fid.variables['yu_ocean'][:]
lon1 = nc_fid.variables['xu_ocean'][:]
z = round(nc_fid.variables['st_ocean'][z_idx]*10)/10


#%% re-organise longitudes
out1 = {'1':'re-arrange data so that maps start at 70W'}
lon_less_m70 = lon1 < -70
lon_less_m70_flipped = np.flipud(lon_less_m70)

for v in var:
    for d, d_out in zip(DATA, DATA_out):
        for m in MMM:
            ma = np.ma.empty((989,3600))
            ma[:,lon_less_m70_flipped] = out[v + d_out + m][:,lon_less_m70]
            ma[:,~lon_less_m70_flipped] = out[v + d_out + m][:,~lon_less_m70]
            out1[v + d_out + m] = ma
            print(v + d_out + m + ' OK !')


#%% prepare data for plot
outf = {'f':'create anomalies for columns 1 and 2'}

for v in var:
    for d, d_out in zip(DATA, DATA_out):
        for m in MMM:
            if d_out is 'CT':
                outf[v + d_out + m] = out1[v + d_out + m]
                
            else:
                outf[v + d_out + m] = out1[v + d_out + m] - out1[v + 'CT' + m]
            
            print(v + d_out + m + ' OK !')
            
            
#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
Bm = Basemap(projection='mill', llcrnrlat=-50.01,urcrnrlat=-24.99,\
llcrnrlon=109.99,urcrnrlon=157.01, resolution='c')

pickle.dump(Bm,open('map.pickle','wb'),-1)  # pickle it 

Bm.fix_aspect = True
m_aspect = Bm.aspect

lon = lon1 + 210

# meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
x, y = Bm(lons, lats)

#
points = np.meshgrid(yy, xx)

row = 3
col = 4

merid = np.zeros(12)
paral = np.zeros(12)

merid[-col:] = 1
paral[np.arange(0,11,col)] = 1

        
#%% plot surface maps: TEMPERATURE
fign = 1
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 14}) 

fig.set_size_inches(11.5, 10) # set figure size in inches
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

        ax = fig.add_subplot(row,col,s)
        pos = ax.get_position()
        bnd = list(pos.bounds)
        magn = 0.04
        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
        ax.set_position(bnd)
        
#        if d_out is 'CT':
        cmap = plt.get_cmap('seismic')
        step = 0.05
        # levels to show on colourbar
        contf_lvls = np.arange(-0.3,0.3+1e-08,step)
        
        #
        ax.set_facecolor('grey')
        pickle.load(open('map.pickle','rb'))   # load here the above pickle
        
        # draw land outlines
        Bm.drawcoastlines(linewidth=0.05)
        Bm.fillcontinents(color='white')
        
        # filled contour plot
        contf = Bm.contourf(x, y, out1['u' + d_out + m], 
                           contf_lvls, cmap=cmap, extend='both')
        
        # for quiver. data points
        yy = np.arange(0, y.shape[0], 20)
        xx = np.arange(0, x.shape[1], 15)
        quiv_scale = 1.5

        # quiver
        Q = Bm.quiver(x[points], y[points],
                 out1['u' + d_out + m][points], 
                 out1['v' + d_out + m][points], 
                 scale=quiv_scale,headwidth=6,
                 headlength=9,headaxislength=9)
        
        # title ...
        if s <= 3:
            title_name = exps[s-1]
            ax.set_title(title_name)
            
        if s is 1 or s is 4 or s is 7 or s is 10:
            ex = ex + 1
            ax.set_ylabel(seasons[ex-1])
            ax.yaxis.labelpad = 35    
        
        # meridians. last input is meridians tick label
        meridians = map(round_to_base, np.arange(110, 180, 10))
        Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,merid[s-1]])
        # parallels. last input is paralles tick label
        parallels = map(round_to_base, np.arange(-45, -10, 5))
        Bm.drawparallels(parallels, linewidth=0.2, labels=[paral[s-1],0,0,0])
        
        if s is 3:
            qk = plt.quiverkey(Q,0.8,1.05,0.1,r'0.1 $m\ s^{-1}$',labelpos='E')
            
        elif s is 12:
            cbar_pos = [bnd[0]+bnd[2]+0.01, bnd[1], 0.01, bnd[3]*4] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='vertical', drawedges=True)
            cbar.set_label(r'Absolute U $(m\ s^{-1})$') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        print('u' + d_out + m + ' OK !')
        
plt.suptitle(r"Sea Surface Velocity", y=0.97, fontsize=20)
#plt.suptitle(r"Sea Surface Temperature ${\theta}$", y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)
if not scriptname:
    scriptname = 'test'
if scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig' + str(fign) + '_' + '.png', 
            bbox_inches='tight', dpi=100)

