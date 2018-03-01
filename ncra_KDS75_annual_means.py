'''
This file creates annual means of Paul's wind experiments in KDS75 in
/g/data3/hh5/tmp/cosima/mom01v5/
and places the outputs in
/g/data/e14/erd561/SAMexp/

Earl Duran 
created: 28-Feb-18
e.duran@unsw.edu.au
'''

import os
#import netCDF4 as nc

input_path = '/g/data3/hh5/tmp/cosima/mom01v5/'
input_folder_path = ['KDS75/', 'KDS75_UP/', 'KDS75_PI/']
file_number = list(range(266, 346, 4))

var = ['temp', 'salt', 'u', 'v']
x_axis = ['xt_ocean', 'xt_ocean', 'xu_ocean', 'xu_ocean']
y_axis = ['yt_ocean', 'yt_ocean', 'yu_ocean', 'yu_ocean']

output_path = '/g/data/e14/erd561/SAMexp/'

for f in input_folder_path:
    input_path_now = input_path + f
    output_path_now = output_path + f

    for n in file_number:

        input_data_path = ''
        for t in range(4):
            input_data_path += ' ' + input_path_now + 'output' + str(n+t) + '/ocean.nc'

        #nc_data = nc.Dataset(input_path_now + 'output' + str(n) + '/ocean.nc', 'w', format='NETCDF4')
        #print(nc_data)
        #time_start = nc_data["temp"]
        time_start = str(n) + '-' + str(n+3)

        for v,x,y in zip(var,x_axis,y_axis):
            #output_data = v + '_' + str(time_start[:]) + '.nc'
            output_data = v + '_' + time_start + '.nc'
            output_data_path = output_path_now + output_data
            pls = os.listdir(output_path_now)

            if output_data not in pls:
                os.system('ncra -d ' + y + ',-60.0,-20.0 -d ' + x + ',-260.0,-190.0 -d st_ocean,0,0 -v ' \
                    + v + input_data_path + ' ' + output_data_path)
                print(output_data_path + ' OK')


