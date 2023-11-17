#!/bin/bash

main=main.cpp
class=class.cu
io=../../src/io/io.cpp

nvcc $io $class $main --std=c++11 -lm -O3 -o acoustic.exe

./acoustic.exe parameters.txt

