#!/bin/bash



ln -sf WFN_parabands.h5 WFN_co.h5
ln -sf WFN_parabands.h5 WFNq_co.h5
srun --ntasks=16   --gpus-per-task=1 kernel.cplx.x &> kernel.inp.out
