#!/bin/bash


mpirun -n 1 ph.x     < dfpt.in &> dfpt.in.out  

python3 ./create_save.py
