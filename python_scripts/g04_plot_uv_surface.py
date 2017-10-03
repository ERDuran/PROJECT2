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
 
MMM = ['JFM', 'AMJ', 'JAS', 'OND']

DATA = ['CT', 'UP', 'SH', 'PI']

var = ['u', 'v']

scriptname = os.path.basename(sys.argv[0])[:-3]

uv0 = {'0':'raw surface u and v data'}

print('OK, all set')


#%% Load data tau_x
z_idx = 0

for v in var:
    for d in DATA:
        for m in MMM:
            # create an instance of the ncCDF4 class
            nc_fid = nc.Dataset(data_path + 'GFDL50_seasons/' + d + '/' + v + 
                                '_' + m + '.nc', 'r')

            # same for temperature (1,2,3)
            uv0[v + d + m] = nc_fid.variables[v][0,z_idx,:,:]
            
            #
            print(v + d + m + ' OK !')
            
# get dimensions
lat = nc_fid.variables['yu_ocean'][:]
lon1 = nc_fid.variables['xu_ocean'][:]
z = nc_fid.variables['st_ocean'][z_idx]


#%% re-organise longitudes
uv1 = {'1':'re-arrange data so that maps start at 70W'}
lon_less_m70 = lon1 < -70 
lon_less_m70_flipped = np.flipud(lon_less_m70)

for v in var:
    for d in DATA:
        for m in MMM:
            ma = np.ma.empty((400,1440))
            ma[:,lon_less_m70_flipped] = uv0[v + d + m][:,lon_less_m70]
            ma[:,~lon_less_m70_flipped] = uv0[v + d + m][:,~lon_less_m70]
            uv1[v + d + m] = ma
            print(v + d + m + ' OK !')


#%% prepare data for plot
uvf = {'f':'create anomalies for columns 1 and 2'}

for v in var:
    for d in DATA:
        for m in MMM:
            if d is 'CT':
                uvf[v + d + m] = uv1[v + d + m]
                
            else:
                uvf[v + d + m] = uv1[v + d + m] - uv1[v + 'CT' + m]
            
            print(v + d + m + ' OK !')
            
            
#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
Bm = Basemap(projection='mill', llcrnrlat=-70,urcrnrlat=0.01,\
llcrnrlon=-70,urcrnrlon=290, resolution='c')

pickle.dump(Bm,open('map.pickle','wb'),-1)  # pickle it 

Bm.fix_aspect = True
m_aspect = Bm.aspect

lon = lon1 + 210


#%% plot surface maps: TEMPERATURE
fign = 1
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 6}) 

fig.set_size_inches(14, 5) # set figure size in inches
row = 4
col = 4

merid = np.zeros(16)
paral = np.zeros(16)

merid[-col:] = 1
paral[np.arange(0,15,col)] = 1

s = 0

for m in MMM:
    for d in DATA:    
        s = s + 1

        # position figure wrt window template
        ax = fig.add_subplot(row,col,s)
        pos = ax.get_position()
        bnd = list(pos.bounds)
        magn = 0.03
        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
        ax.set_position(bnd)
        
        # colourmap: blue to red because looking at bias.
        if d is 'CT':
            cmap = plt.get_cmap('RdYlGn_r')
            step = 0.05
            # levels to show on colourbar
            contf_lvls = np.arange(-0.25,0.25+1e-08,step)
            
        else:
            cmap = plt.get_cmap('seismic')
            step = 0.02
            # levels to show on colourbar
            contf_lvls = np.arange(-0.14,0.14+1e-08,step)            
        
        
        ax.set_facecolor('grey')
        
        # meshgrid of lon lats
        lons, lats = np.meshgrid(lon, lat)
        # use projection template to create x and y axis
        x, y = Bm(lons, lats)
        
        # for quiver. data points
        yy = np.arange(0, y.shape[0], 20)
        xx = np.arange(0, x.shape[1], 40)
        quiv_scale = 5
        
        #
        points = np.meshgrid(yy, xx)
        
        #
        pickle.load(open('map.pickle','rb'))   # load here the above pickle
        
        # draw land outlines
        Bm.drawcoastlines(linewidth=0.05)
        Bm.fillcontinents(color='white')
        
        # filled contour plot
        contf = Bm.contourf(x, y, uvf['u' + d + m], 
                           contf_lvls, cmap=cmap, extend='both')
        
        if d is 'CT':
            # quiver
            Bm.quiver(x[points], y[points],
                      uvf['u' + d + m][points], 
                      uvf['v' + d + m][points], 
                      scale=quiv_scale)
        
        # title ...
        title_name = 'GFDL50 ' + d + ' u ' + m + ' z = ' + str(z)
        ax.set_title(title_name)
        
        def round_to_base(x, base=5):
            return int(base * round(float(x) / base))
        # meridians. last input is meridians tick label
        meridians = map(round_to_base, np.arange(-160, 180, 40))
        Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,merid[s-1]])
        # parallels. last input is paralles tick label
        parallels = map(round_to_base, np.arange(-80, 20, 20))
        Bm.drawparallels(parallels, linewidth=0.2, labels=[paral[s-1],0,0,0])
        
        if s is 13:
            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='horizontal', drawedges=True)
            cbar.set_label('m/s') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        elif s is 14:
            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='horizontal', drawedges=True)
            cbar.set_label('m/s') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        print(d + m + ' OK !')
        

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig' + str(fign) + '_' + str(z_idx) + '.png',
            bbox_inches='tight', dpi=300)


#%% plot surface maps: TEMPERATURE
fign = 2
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 6}) 

fig.set_size_inches(14, 5) # set figure size in inches
row = 4
col = 4

merid = np.zeros(16)
paral = np.zeros(16)

merid[-col:] = 1
paral[np.arange(0,15,col)] = 1

s = 0

for m in MMM:
    for d in DATA:    
        s = s + 1

        # position figure wrt window template
        ax = fig.add_subplot(row,col,s)
        pos = ax.get_position()
        bnd = list(pos.bounds)
        magn = 0.03
        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
        ax.set_position(bnd)
        
        # colourmap: blue to red because looking at bias.
        cmap = plt.get_cmap('RdYlGn_r')
        step = 0.05
        # levels to show on colourbar
        contf_lvls = np.arange(-0.25,0.25+1e-08,step)
        
        
        ax.set_facecolor('grey')
        
        # meshgrid of lon lats
        lons, lats = np.meshgrid(lon, lat)
        # use projection template to create x and y axis
        x, y = Bm(lons, lats)
        
        # for quiver. data points
        yy = np.arange(0, y.shape[0], 20)
        xx = np.arange(0, x.shape[1], 40)
        quiv_scale = 5
        
        #
        points = np.meshgrid(yy, xx)
        
        #
        pickle.load(open('map.pickle','rb'))   # load here the above pickle
        
        # draw land outlines
        Bm.drawcoastlines(linewidth=0.05)
        Bm.fillcontinents(color='white')
        
        # filled contour plot
        contf = Bm.contourf(x, y, uv1['u' + d + m], 
                           contf_lvls, cmap=cmap, extend='both')
        
        # quiver
        Bm.quiver(x[points], y[points],
                  uv1['u' + d + m][points], 
                  uv1['v' + d + m][points], 
                  scale=quiv_scale)
        
        # title ...
        title_name = 'GFDL50 ' + d + ' u ' + m + ' z = ' + str(z)
        ax.set_title(title_name)
        
        def round_to_base(x, base=5):
            return int(base * round(float(x) / base))
        # meridians. last input is meridians tick label
        meridians = map(round_to_base, np.arange(-160, 180, 40))
        Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,merid[s-1]])
        # parallels. last input is paralles tick label
        parallels = map(round_to_base, np.arange(-80, 20, 20))
        Bm.drawparallels(parallels, linewidth=0.2, labels=[paral[s-1],0,0,0])
        
        if s is 13:
            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='horizontal', drawedges=True)
            cbar.set_label('m/s') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        print(d + m + ' OK !')
        

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig' + str(fign) + '_' + str(z_idx) + '.png',
            bbox_inches='tight', dpi=300)


#%% plot surface maps: TEMPERATURE
fign = 3
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 6}) 

fig.set_size_inches(14, 5) # set figure size in inches
row = 4
col = 4

merid = np.zeros(16)
paral = np.zeros(16)

merid[-col:] = 1
paral[np.arange(0,15,col)] = 1

s = 0

for m in MMM:
    for d in DATA:    
        s = s + 1

        # position figure wrt window template
        ax = fig.add_subplot(row,col,s)
        pos = ax.get_position()
        bnd = list(pos.bounds)
        magn = 0.03
        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
        ax.set_position(bnd)
        
        # colourmap: blue to red because looking at bias.
        if d is 'CT':
            cmap = plt.get_cmap('RdYlGn_r')
            step = 0.05
            # levels to show on colourbar
            contf_lvls = np.arange(-0.25,0.25+1e-08,step)
            
        else:
            cmap = plt.get_cmap('seismic')
            step = 0.02
            # levels to show on colourbar
            contf_lvls = np.arange(-0.14,0.14+1e-08,step)            
        
        
        ax.set_facecolor('grey')
        
        # meshgrid of lon lats
        lons, lats = np.meshgrid(lon, lat)
        # use projection template to create x and y axis
        x, y = Bm(lons, lats)
        
        # for quiver. data points
        yy = np.arange(0, y.shape[0], 20)
        xx = np.arange(0, x.shape[1], 40)
        quiv_scale = 5
        
        #
        points = np.meshgrid(yy, xx)
        
        #
        pickle.load(open('map.pickle','rb'))   # load here the above pickle
        
        # draw land outlines
        Bm.drawcoastlines(linewidth=0.05)
        Bm.fillcontinents(color='white')
        
        # filled contour plot
        contf = Bm.contourf(x, y, uvf['v' + d + m], 
                           contf_lvls, cmap=cmap, extend='both')
        
        if d is 'CT':
            # quiver
            Bm.quiver(x[points], y[points],
                      uvf['u' + d + m][points], 
                      uvf['v' + d + m][points], 
                      scale=quiv_scale)
        
        # title ...
        title_name = 'GFDL50 ' + d + ' v ' + m + ' z = ' + str(z)
        ax.set_title(title_name)
        
        def round_to_base(x, base=5):
            return int(base * round(float(x) / base))
        # meridians. last input is meridians tick label
        meridians = map(round_to_base, np.arange(-160, 180, 40))
        Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,merid[s-1]])
        # parallels. last input is paralles tick label
        parallels = map(round_to_base, np.arange(-80, 20, 20))
        Bm.drawparallels(parallels, linewidth=0.2, labels=[paral[s-1],0,0,0])
        
        if s is 13:
            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='horizontal', drawedges=True)
            cbar.set_label('m/s') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        elif s is 14:
            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='horizontal', drawedges=True)
            cbar.set_label('m/s') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        print(d + m + ' OK !')
        

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig' + str(fign) + '_' + str(z_idx) + '.png',
            bbox_inches='tight', dpi=300)


#%% plot surface maps: TEMPERATURE
fign = 4
plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
matplotlib.rcParams.update({'font.size': 6}) 

fig.set_size_inches(14, 5) # set figure size in inches
row = 4
col = 4

merid = np.zeros(16)
paral = np.zeros(16)

merid[-col:] = 1
paral[np.arange(0,15,col)] = 1

s = 0

for m in MMM:
    for d in DATA:    
        s = s + 1

        # position figure wrt window template
        ax = fig.add_subplot(row,col,s)
        pos = ax.get_position()
        bnd = list(pos.bounds)
        magn = 0.03
        bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
        ax.set_position(bnd)
        
        # colourmap: blue to red because looking at bias.
        cmap = plt.get_cmap('RdYlGn_r')
        step = 0.05
        # levels to show on colourbar
        contf_lvls = np.arange(-0.25,0.25+1e-08,step)
        
        
        ax.set_facecolor('grey')
        
        # meshgrid of lon lats
        lons, lats = np.meshgrid(lon, lat)
        # use projection template to create x and y axis
        x, y = Bm(lons, lats)
        
        # for quiver. data points
        yy = np.arange(0, y.shape[0], 20)
        xx = np.arange(0, x.shape[1], 40)
        quiv_scale = 5
        
        #
        points = np.meshgrid(yy, xx)
        
        #
        pickle.load(open('map.pickle','rb'))   # load here the above pickle
        
        # draw land outlines
        Bm.drawcoastlines(linewidth=0.05)
        Bm.fillcontinents(color='white')
        
        # filled contour plot
        contf = Bm.contourf(x, y, uv1['v' + d + m], 
                           contf_lvls, cmap=cmap, extend='both')
        
        # quiver
        Bm.quiver(x[points], y[points],
                  uv1['u' + d + m][points], 
                  uv1['v' + d + m][points], 
                  scale=quiv_scale)
        
        # title ...
        title_name = 'GFDL50 ' + d + ' v ' + m + ' z = ' + str(z)
        ax.set_title(title_name)
        
        def round_to_base(x, base=5):
            return int(base * round(float(x) / base))
        # meridians. last input is meridians tick label
        meridians = map(round_to_base, np.arange(-160, 180, 40))
        Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,merid[s-1]])
        # parallels. last input is paralles tick label
        parallels = map(round_to_base, np.arange(-80, 20, 20))
        Bm.drawparallels(parallels, linewidth=0.2, labels=[paral[s-1],0,0,0])
        
        if s is 13:
            cbar_pos = [bnd[0], bnd[1]-0.03, bnd[2], 0.01] 
            cbar_axe = fig.add_axes(cbar_pos)
            cbar = plt.colorbar(contf, cax=cbar_axe, 
                                orientation='horizontal', drawedges=True)
            cbar.set_label('m/s') # units label on colourbar
            cbar.dividers.set_linewidth(0.2)
            cbar.outline.set_linewidth(0.5)
            cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
            cbar.ax.tick_params(width=0.2, length= 2)
            
        print(d + m + ' OK !')
        

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig' + str(fign) + '_' + str(z_idx) + '.png',
            bbox_inches='tight', dpi=300)
