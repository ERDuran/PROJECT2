#!/bin/bash

input_path_orig=GFDL50/

ptdm_files=553

input_path_next=(CT SH UP PI)

var=(tau_x tau_y temp salt u v)

for d in {0..3}; do
    input_path=${input_path_orig}${input_path_next[${d}]}/

    for v in {0..5}; do
		for y in {0..7}; do
		    JFM_data[${y}]=${input_path}${var[${v}]}_JFM_$(expr ${ptdm_files} + ${y}).nc
		    AMJ_data[${y}]=${input_path}${var[${v}]}_AMJ_$(expr ${ptdm_files} + ${y}).nc
		    JAS_data[${y}]=${input_path}${var[${v}]}_JAS_$(expr ${ptdm_files} + ${y}).nc
		    OND_data[${y}]=${input_path}${var[${v}]}_OND_$(expr ${ptdm_files} + ${y}).nc
		done

		ncra -O ${JFM_data[@]} GFDL50_seasons/${input_path_next[${d}]}/${var[${v}]}_JFM.nc
		ncra -O ${AMJ_data[@]} GFDL50_seasons/${input_path_next[${d}]}/${var[${v}]}_AMJ.nc
		ncra -O ${JAS_data[@]} GFDL50_seasons/${input_path_next[${d}]}/${var[${v}]}_JAS.nc
		ncra -O ${OND_data[@]} GFDL50_seasons/${input_path_next[${d}]}/${var[${v}]}_OND.nc

		echo "${input_path_next[${d}]} ${var[${v}]} OK !"
    done
done