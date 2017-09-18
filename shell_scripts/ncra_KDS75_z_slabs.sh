#!/bin/bash

input_path_orig=/g/data3/hh5/tmp/cosima/mom01v5/

MMM_files=(314 318 322 326 330 334 338 342)

input_path_next=(KDS75 KDS75_UP KDS75_PI)

var1=(tau_x tau_y)

var2=(temp salt u v)

y_axis=(yt_ocean yt_ocean yu_ocean yu_ocean)

z_min=(0 15 30 45)
z_max=(14 29 44 59)

for d in {0..2}; do
    input_path=${input_path_orig}${input_path_next[${d}]}/

    for t in {0..7}; do
        JFM_data1[${t}]=${input_path}output${MMM_files[${t}]}/ocean_month.nc
	    AMJ_data1[${t}]=${input_path}output$(expr ${MMM_files[${t}]} + 1)/ocean_month.nc
        JAS_data1[${t}]=${input_path}output$(expr ${MMM_files[${t}]} + 2)/ocean_month.nc
        OND_data1[${t}]=${input_path}output$(expr ${MMM_files[${t}]} + 3)/ocean_month.nc

        JFM_data2[${t}]=${input_path}output${MMM_files[${t}]}/ocean.nc
        AMJ_data2[${t}]=${input_path}output$(expr ${MMM_files[${t}]} + 1)/ocean.nc
        JAS_data2[${t}]=${input_path}output$(expr ${MMM_files[${t}]} + 2)/ocean.nc
        OND_data2[${t}]=${input_path}output$(expr ${MMM_files[${t}]} + 3)/ocean.nc
    done

    for v1 in {0..1}; do
    	ncra -O -d time,0,0 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${JFM_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Jan.nc
    	ncra -O -d time,1,1 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${JFM_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Feb.nc
    	ncra -O -d time,2,2 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${JFM_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Mar.nc
    	ncra -O -d time,0,0 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${AMJ_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Apr.nc
    	ncra -O -d time,1,1 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${AMJ_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_May.nc
    	ncra -O -d time,2,2 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${AMJ_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Jun.nc
    	ncra -O -d time,0,0 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${JAS_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Jul.nc
    	ncra -O -d time,1,1 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${JAS_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Aug.nc
    	ncra -O -d time,2,2 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${JAS_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Sep.nc
    	ncra -O -d time,0,0 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${OND_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Oct.nc
    	ncra -O -d time,1,1 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${OND_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Nov.nc
    	ncra -O -d time,2,2 -d yu_ocean,-70.2,0.2 -v ${var1[${v1}]} ${OND_data1[@]} ${input_path_next[${d}]}/${var1[${v1}]}_Dec.nc

    	echo "${input_path_next[${d}]} ${var1[${v1}]} OK"
    done

    for v2 in {0..3}; do
    	for z in {0..3}; do
            ncra -O -d ${y_axis[${v2}]},-70.2,0.2 -d st_ocean,${z_min[${z}]},${z_max[${z}]} -v ${var2[${v2}]} ${JFM_data2[@]} ${input_path_next[${d}]}/${var2[${v2}]}_JFM_${z}.nc
            echo "JFM OK"
            ncra -O -d ${y_axis[${v2}]},-70.2,0.2 -d st_ocean,${z_min[${z}]},${z_max[${z}]} -v ${var2[${v2}]} ${AMJ_data2[@]} ${input_path_next[${d}]}/${var2[${v2}]}_AMJ_${z}.nc
            echo "AMJ OK"
            ncra -O -d ${y_axis[${v2}]},-70.2,0.2 -d st_ocean,${z_min[${z}]},${z_max[${z}]} -v ${var2[${v2}]} ${JAS_data2[@]} ${input_path_next[${d}]}/${var2[${v2}]}_JAS_${z}.nc
            echo "JAS OK"
            ncra -O -d ${y_axis[${v2}]},-70.2,0.2 -d st_ocean,${z_min[${z}]},${z_max[${z}]} -v ${var2[${v2}]} ${OND_data2[@]} ${input_path_next[${d}]}/${var2[${v2}]}_OND_${z}.nc
            echo "OND OK"

            echo "${input_path_next[${d}]} ${var2[${v2}]} z=${z} OK"
        done
    done
done

