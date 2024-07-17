#!/bin/bash


mpirun -n 1 projwfc.x -pd .true. < pdos.in &> pdos.in.out 
