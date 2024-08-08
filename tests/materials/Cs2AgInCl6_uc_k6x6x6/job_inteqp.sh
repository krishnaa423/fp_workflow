#!/bin/bash



ln -sf WFN_coo ./WFN_co 
ln -sf WFN_dftelbands ./WFN_fi 
ln -sf ./eqp1.dat ./eqp_co.dat 
srun --ntasks=16   --gpus-per-task=1 inteqp.cplx.x &> inteqp.inp.out 
mv bandstructure.dat bandstructure_inteqp.dat 
