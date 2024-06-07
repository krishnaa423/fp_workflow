#!/bin/bash 

cd ../../utils
python3 conf_bse.py --run &> ../convergence/dft/conv_dft.out 
cd ../convergence/dft/