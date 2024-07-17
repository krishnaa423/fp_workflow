#!/bin/bash


mpirun -n 1 pw.x < dftelbands.in &> dftelbands.in.out  
cp ./tmp/struct.xml ./dftelbands.xml 
