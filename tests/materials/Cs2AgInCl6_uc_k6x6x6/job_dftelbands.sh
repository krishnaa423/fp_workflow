#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 pw.x < dftelbands.in &> dftelbands.in.out  
cp ./tmp/struct.xml ./dftelbands.xml 
