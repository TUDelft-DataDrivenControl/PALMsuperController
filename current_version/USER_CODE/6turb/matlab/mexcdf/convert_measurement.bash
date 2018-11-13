#!/bin/bash
Wpname="6turb"
nturbines=6

if [ -f "JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters01" ]; then
for i in $(seq 1 $nturbines ); do
	# store the turbine file as .txt file
    iconv --from-code US-ASCII --to-code UTF-8 -c JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i} > JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}.txt # not suitable for 10+ turbines yet
 	# remove first two line   
    sed -i 1,2d JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}.txt
    # remove original file
    rm JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}
done
fi

rm JOBS/$Wpname/MONITORING/${Wpname}_cpu
rm JOBS/$Wpname/MONITORING/${Wpname}_header
rm JOBS/$Wpname/MONITORING/${Wpname}_rc
