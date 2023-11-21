#!/bin/bash

io=../../src/io/io.cpp
main=cppCodes/scalar_wave_main.cpp
class=cppCodes/scalar_wave_class.cu

cd ../../run

./program.sh -compile; clear 
./program.sh -modeling; clear

cd ../tools/seismogram

nvcc $io $class $main --std=c++11 -lm -O3 -o acoustic.exe

./acoustic.exe parameters.txt; clear  

rm acoustic.exe

python3 pyCodes/add_trace_header.py

rm segy_data/*.bin
