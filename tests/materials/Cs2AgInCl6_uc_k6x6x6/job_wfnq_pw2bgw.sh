#!/bin/bash



srun --ntasks=4   --gpus-per-task=1 pw2bgw.x -pd .true. < wfnq_pw2bgw.in &> wfnq_pw2bgw.in.out 
cp ./tmp/WFNq_coo ./
wfn2hdf.x BIN WFNq_coo WFNq_coo.h5 
