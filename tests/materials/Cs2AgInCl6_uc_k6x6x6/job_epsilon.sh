#!/bin/bash



ln -sf WFN_parabands.h5 ./WFN.h5 
ln -sf WFNq_coo.h5 ./WFNq.h5 
srun --ntasks=16   --gpus-per-task=1 epsilon.cplx.x &> epsilon.inp.out 