#!/bin/bash


mpirun -n 1 pw2bgw.x -pd .true. < wfn_pw2bgw.in &> wfn_pw2bgw.in.out 
cp ./tmp/WFN_coo ./
cp ./tmp/RHO ./
cp ./tmp/VXC ./
cp ./tmp/VSC ./
cp ./tmp/VKB ./
