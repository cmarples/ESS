#!/bin/bash
# Run from ESS directory

# Move to build directory and compile
cd build
rm -rf *
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
cmake --build .

# Run the program to generate data
cd ..
./build/run_sampling < input_area.txt
echo "Patch Areas Calculated"
./build/run_sampling < input_time.txt
./build/run_sampling < input_time2.txt
echo "Timings Measured"
./build/run_sampling < input_dist.txt
./build/run_sampling < input_dist2.txt
echo "Distributions Built"
./build/run_sampling < input_samp.txt
echo "Samples Generated"

echo "Program Complete"
