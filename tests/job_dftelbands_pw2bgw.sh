#!/bin/bash


mpirun -n 1 pw2bgw.x -pd .true. < dftelbands_pw2bgw.in &> dftelbands_pw2bgw.in.out 
cp ./tmp/WFN_dftelbands ./
wfn2hdf.x BIN WFN_dftelbands WFN_dftelbands.h5
