#!/bin/bash

declare -a trace=(../data/traj_pbf/region5C-large.trace-ss20-n5.pbf)

declare -a bbox=(../data/bbox/region5C-large.trace-ss20-n5.yaml)

declare -a log=(log_region5C-large.trace-ss20-n5.txt)


 
for i in `seq 0 0`; 
do 	 
	echo ${log[$i]}
	bin/GPStudio --interactive=false --trace_fname=${trace[$i]}  --bbox_fname=${bbox[i]}   2>&1 | tee ${log[$i]}

	source makeOutput.sh >> ${log[$i]}
	echo ${trace[$i]} >> ${log[$i]}
	tail -n 25 ${log[$i]}
done


exit

