#!/bin/bash

echo Enter wind farm:
read Wpname

echo Enter number of turbines:
read nturbines

echo Restart job y or n:
read pre

if [ -f "JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters01" ]; then
for i in $(seq 1 $nturbines ); do
	
	if [ "$i" -lt 10 ]; then
	# store the turbine file as .txt file
    iconv --from-code US-ASCII --to-code UTF-8 -c JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i} > JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}.txt 
 	# remove first two line   
    sed -i 1,2d JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}.txt
    # move turbine files
    mv JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}.txt COPIED_FILES/
    # remove original file
    rm JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters0${i}

    else
	# store the turbine file as .txt file
    iconv --from-code US-ASCII --to-code UTF-8 -c JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters${i} > JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters${i}.txt 
 	# remove first two line   
    sed -i 1,2d JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters${i}.txt
    # move turbine files
    mv JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters${i}.txt COPIED_FILES/
    # remove original file
    rm JOBS/$Wpname/MONITORING/${Wpname}_turbine_parameters${i}    
   
    fi
	 
done
fi

rm JOBS/$Wpname/MONITORING/${Wpname}_cpu
rm JOBS/$Wpname/MONITORING/${Wpname}_header
rm JOBS/$Wpname/MONITORING/${Wpname}_rc
mv JOBS/$Wpname/OUTPUT/${Wpname}_m01.nc COPIED_FILES/
mv ../../job_queue/lchpc06_${Wpname}* COPIED_FILES/ 
rm ../../job_queue/${Wpname}* 

# for restart run
if [ "$pre" = "y" ]; then
 cp JOBS/$Wpname/INPUT/${Wpname}_p3df COPIED_FILES/
else
 cp JOBS/$Wpname/INPUT/${Wpname}_p3d COPIED_FILES/	 
fi
