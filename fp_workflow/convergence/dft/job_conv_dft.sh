#!/bin/bash 

cd ../../utils
python3 conf_dft.py --run &> ../convergence/dft/conv_dft.out 
cd ../convergence/dft/