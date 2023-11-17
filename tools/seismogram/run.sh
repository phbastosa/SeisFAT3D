#!/bin/bash

io=../../src/io/io.cpp
main=scalar_wave_main.cpp
class=scalar_wave_class.cu

nvcc $io $class $main --std=c++11 -lm -O3 -o acoustic.exe

./acoustic.exe parameters.txt

rm acoustic.exe  

python3 add_trace_header.py

rm segy_data/*.bin

