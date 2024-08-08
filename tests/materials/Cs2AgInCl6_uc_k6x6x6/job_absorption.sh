#!/bin/bash



ln -sf WFN_parabands.h5 WFN_co.h5 
ln -sf WFN_parabands.h5 WFNq_co.h5 
ln -sf WFN_parabands.h5 WFN_fi.h5 
ln -sf WFN_parabands.h5 WFNq_fi.h5 
ln -sf eqp1.dat eqp_co.dat 
srun --ntasks=16   --gpus-per-task=1 absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
