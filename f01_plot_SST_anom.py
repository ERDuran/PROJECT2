'''
This file plots annual means SST from Paul's wind experiments in KDS75 in
/Users/earl/Desktop/Yang
and places the outputs in
/Users/earl/Dropbox/SAMexp/figures

*** Run Yang.sh to mount data !! ***

Earl Duran 
created: 28-Feb-18
e.duran@unsw.edu.au
'''

import os
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import sys
import pickle
def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()
#def round_to_base(x, base=5):
    #return int(base * round(float(x) / base))

scriptname = os.path.basename(sys.argv[0])[:-3]

data_path = '/Users/earl/Desktop/Yang/'
figures_path = '/Users/earl/Dropbox/SAMexp/figures/'

input_folder_path = ['KDS75/', 'KDS75_UP/', 'KDS75_PI/']
file_number = list(range(266, 346, 4))
var = ['temp', 'salt', 'u', 'v']

CT_temp = {'year':'output'}
PI_temp = {'year':'output'}
PI_temp_anom = {'year':'output'}
time = []

for f in file_number:
    CT = nc.Dataset(data_path + input_folder_path[0] + var[0] + '_' + str(f) + '-' + str(f+3) + '.nc', 'r')
    PI = nc.Dataset(data_path + input_folder_path[2] + var[0] + '_' + str(f) + '-' + str(f+3) + '.nc', 'r')

    time_now = int(round(CT.variables['time'][0]/365))

    time.append(time_now)
    
    CT_temp[str(time_now)] = CT.variables[var[0]][0,0,:,:]
    PI_temp[str(time_now)] = PI.variables[var[0]][0,0,:,:]

    PI_temp_anom[str(time_now)] = PI_temp[str(time_now)] - CT_temp[str(time_now)]
    print(time_now)


lat = CT.variables['yt_ocean'][:]
lon1 = CT.variables['xt_ocean'][:]
lon = lon1 + 360



            
#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
bm = Basemap(projection='mill', llcrnrlat=-60,urcrnrlat=-20,\
llcrnrlon=100,urcrnrlon=170, resolution='c')

bm.fix_aspect = True
m_aspect = bm.aspect

#%% plot surface maps: TEMPERATURE
matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 2

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for t in time:
    plt.close('all') # close all existing figures
    fig = plt.figure() # generate figure
    fig.set_size_inches(16, 8) # set figure size in inches
    
    # position figure wrt window template
    ax = fig.add_subplot(row,col,1)
    pos = ax.get_position()
    bnd = list(pos.bounds)
    magn = 0.04
    bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
    ax.set_position(bnd)
    
    # colourmap: blue to red because looking at bias.
    cmap = plt.get_cmap('gist_ncar')
    step = 2
    # levels to show on colourbar
    contf_lvls = np.arange(-2,30+1e-08,step)
        
    #cmap = plt.get_cmap('seismic')
    #step = 0.02
    ## levels to show on colourbar
    #contf_lvls = np.arange(-0.14,0.14+1e-08,step)            
    
    
    ax.set_facecolor('grey')
    
     # meshgrid of lon lats
    lons, lats = np.meshgrid(lon, lat)
    # use projection template to create x and y axis
    Bm_lons, Bm_lats = bm(lons, lats)   
            
    
    # draw land outlines
    bm.drawcoastlines(linewidth=0.05)
    bm.fillcontinents(color='white')
    
    # filled contour plot
    contf = bm.contourf(Bm_lons, Bm_lats, PI_temp[str(t)],contf_lvls, cmap=cmap, extend='both')
    #contf = bm.contourf(Bm_lons, Bm_lats, PI_temp_anom[str(t)])
    
       
    # title ...
    ax.set_title('PI year ' + str(t))
    
    # meridians. last input is meridians tick label
    bm.drawmeridians(np.arange(110, 180, 10), linewidth=0.2, labels=[0,0,0,1])
    # parallels. last input is paralles tick label
    bm.drawparallels(np.arange(-60, -20, 5), linewidth=0.2, labels=[1,0,0,0])
      
    #cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
    #cbar_axe = fig.add_axes(cbar_pos)
    cbar = plt.colorbar(contf, orientation='horizontal', drawedges=True)
    cbar.set_label(r'SST $^{\circ}C$') # units label on colourbar
    cbar.dividers.set_linewidth(0.2)
    cbar.outline.set_linewidth(0.5)
    cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
    cbar.ax.tick_params(width=0.2, length= 2)
    #cbar.ax.yaxis.set_label_position('left')
    #cbar.ax.yaxis.set_ticks_position('left')
    

    ###########################################################################
    ax = fig.add_subplot(row,col,2)
    pos = ax.get_position()
    bnd = list(pos.bounds)
    magn = 0.04
    bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
    ax.set_position(bnd)
    
    # colourmap: blue to red because looking at bias.
    cmap = plt.get_cmap('seismic')
    step = 0.25
    # levels to show on colourbar
    contf_lvls = np.arange(-2,2+1e-08,step)
        
    #cmap = plt.get_cmap('seismic')
    #step = 0.02
    ## levels to show on colourbar
    #contf_lvls = np.arange(-0.14,0.14+1e-08,step)            
    
    
    ax.set_facecolor('grey')
    
     # meshgrid of lon lats
    lons, lats = np.meshgrid(lon, lat)
    # use projection template to create x and y axis
    Bm_lons, Bm_lats = bm(lons, lats)   
            
    
    # draw land outlines
    bm.drawcoastlines(linewidth=0.05)
    bm.fillcontinents(color='white')
    
    # filled contour plot
    contf = bm.contourf(Bm_lons, Bm_lats, PI_temp_anom[str(t)],contf_lvls, cmap=cmap, extend='both')
    #contf = bm.contourf(Bm_lons, Bm_lats, PI_temp_anom[str(t)])
    
       
    # title ...
    ax.set_title('PI - CT year ' + str(t) + '. CT is concomitantly extended')
    
    # meridians. last input is meridians tick label
    bm.drawmeridians(np.arange(110, 180, 10), linewidth=0.2, labels=[0,0,0,1])
    # parallels. last input is paralles tick label
    bm.drawparallels(np.arange(-60, -20, 5), linewidth=0.2, labels=[1,0,0,0])
      
    #cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
    #cbar_axe = fig.add_axes(cbar_pos)
    cbar = plt.colorbar(contf, orientation='horizontal', drawedges=True)
    cbar.set_label(r'SST anomaly $^{\circ}C$') # units label on colourbar
    cbar.dividers.set_linewidth(0.2)
    cbar.outline.set_linewidth(0.5)
    cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
    cbar.ax.tick_params(width=0.2, length= 2)
    #cbar.ax.yaxis.set_label_position('left')
    #cbar.ax.yaxis.set_ticks_position('left')

        
    print(str(t) + ' OK !')
            
    plt.suptitle(r"KDS75 year to year SST in PI which starts at year 69.", 
                 y=0.97, fontsize=20)
    
    output_ls = os.listdir(figures_path)
    
    #
    if not scriptname:
        scriptname = 'test'
    
    elif scriptname not in output_ls:
        os.mkdir(figures_path + '/' + scriptname)
    
    
    # save figure
    plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
                + '_fig1_' + str(t) + '.png', bbox_inches='tight', dpi=200)

