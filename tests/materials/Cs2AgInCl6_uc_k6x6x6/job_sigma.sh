#!/bin/bash



ln -sf WFN_parabands.h5 ./WFN_inner.h5 
srun --ntasks=16   --gpus-per-task=1 sigma.cplx.x &> sigma.inp.out
