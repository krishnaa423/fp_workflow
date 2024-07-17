#!/bin/bash


ln -sf WFN_parabands.h5 ./WFN_inner.h5 
mpirun -n 1 sigma.cplx.x &> sigma.inp.out
