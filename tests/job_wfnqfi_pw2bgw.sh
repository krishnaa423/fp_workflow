#!/bin/bash


mpirun -n 1 pw2bgw.x -pd .true. < wfnqfi_pw2bgw.in &> wfnqfi_pw2bgw.in.out 
cp ./tmp/WFNq_fii ./
wfn2hdf.x BIN WFNq_fii WFNq_fii.h5 
