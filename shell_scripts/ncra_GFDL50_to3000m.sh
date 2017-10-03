#!/bin/bash

input_path_orig=/g/data3/hh5/tmp/pas561/

ptdm_files=553

JFM=($(seq -w 003 5 088))
AMJ=($(seq -w 093 5 178))
JAS=($(seq -w 183 5 268))
OND=($(seq -w 273 5 358))

input_path_next=(gfdl_nyf_1080_cp gfdl_nyf_1080_SH_ICyr200 gfdl_nyf_1080_UP gfdl_nyf_1080_PI_ICyr200)

output_path=(CT SH UP PI)

var1=(tau_x tau_y)

var2=(temp salt u v)

y_axis=(yt_ocean yt_ocean yu_ocean yu_ocean)

for d in {0..3}; do
    input_path=${input_path_orig}${input_path_next[${d}]}/

    for y in {0..7}; do
        for t in {0..17}; do
            JFM_data[${t}]=${input_path}output$(expr ${ptdm_files} + ${y})/ocean__$(expr ${ptdm_files} + ${y})_${JFM[${t}]}.nc
    	    AMJ_data[${t}]=${input_path}output$(expr ${ptdm_files} + ${y})/ocean__$(expr ${ptdm_files} + ${y})_${AMJ[${t}]}.nc
            JAS_data[${t}]=${input_path}output$(expr ${ptdm_files} + ${y})/ocean__$(expr ${ptdm_files} + ${y})_${JAS[${t}]}.nc
            OND_data[${t}]=${input_path}output$(expr ${ptdm_files} + ${y})/ocean__$(expr ${ptdm_files} + ${y})_${OND[${t}]}.nc
        done

        for v1 in {0..1}; do
        	ncra -O -d yu_ocean,-70.5,0.5 -v ${var1[${v1}]} ${JFM_data[@]} GFDL50/${output_path[${d}]}/${var1[${v1}]}_JFM_$(expr ${ptdm_files} + ${y}).nc
        	ncra -O -d yu_ocean,-70.5,0.5 -v ${var1[${v1}]} ${AMJ_data[@]} GFDL50/${output_path[${d}]}/${var1[${v1}]}_AMJ_$(expr ${ptdm_files} + ${y}).nc
        	ncra -O -d yu_ocean,-70.5,0.5 -v ${var1[${v1}]} ${JAS_data[@]} GFDL50/${output_path[${d}]}/${var1[${v1}]}_JAS_$(expr ${ptdm_files} + ${y}).nc
        	ncra -O -d yu_ocean,-70.5,0.5 -v ${var1[${v1}]} ${OND_data[@]} GFDL50/${output_path[${d}]}/${var1[${v1}]}_OND_$(expr ${ptdm_files} + ${y}).nc

        	echo "${output_path[${d}]} ${var1[${v1}]} $(expr ${ptdm_files} + ${y}) OK !"
        done

        for v2 in {0..3}; do
            ncra -O -d ${y_axis[${v2}]},-70.5,0.5 -d st_ocean,0,37 -v ${var2[${v2}]} ${JFM_data[@]} GFDL50/${output_path[${d}]}/${var2[${v2}]}_JFM_$(expr ${ptdm_files} + ${y}).nc
            echo "JFM OK"
            ncra -O -d ${y_axis[${v2}]},-70.5,0.5 -d st_ocean,0,37 -v ${var2[${v2}]} ${AMJ_data[@]} GFDL50/${output_path[${d}]}/${var2[${v2}]}_AMJ_$(expr ${ptdm_files} + ${y}).nc
            echo "AMJ OK"
            ncra -O -d ${y_axis[${v2}]},-70.5,0.5 -d st_ocean,0,37 -v ${var2[${v2}]} ${JAS_data[@]} GFDL50/${output_path[${d}]}/${var2[${v2}]}_JAS_$(expr ${ptdm_files} + ${y}).nc
            echo "JAS OK"
            ncra -O -d ${y_axis[${v2}]},-70.5,0.5 -d st_ocean,0,37 -v ${var2[${v2}]} ${OND_data[@]} GFDL50/${output_path[${d}]}/${var2[${v2}]}_OND_$(expr ${ptdm_files} + ${y}).nc
            echo "OND OK"

            echo "${output_path[${d}]} ${var2[${v2}]} $(expr ${ptdm_files} + ${y}) OK !"
        done
    done
done