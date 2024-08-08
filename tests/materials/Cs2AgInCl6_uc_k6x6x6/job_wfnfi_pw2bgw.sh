#!/bin/bash



srun --ntasks=4   --gpus-per-task=1 pw2bgw.x -pd .true. < wfnfi_pw2bgw.in &> wfnfi_pw2bgw.in.out 
cp ./tmp/WFN_fii ./
wfn2hdf.x BIN WFN_fii WFN_fii.h5 
