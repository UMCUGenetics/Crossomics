#!/bin/bash
# This script is to run the run_array.sh script with all seeds


for seed in 2528 6140 3880 2771 8455 3200 6250 4860 6297 244 3764 2464 3218 2282 5600 2359 8353 6399 2001
do
	echo "queued $seed"
	./run_array.sh $seed
done
