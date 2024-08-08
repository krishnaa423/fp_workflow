#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 parabands.cplx.x &> parabands.inp.out 
