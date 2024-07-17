#!/bin/bash


mpirun -n 1 projwfc.x -pd .true. < kpdos.in &> kpdos.in.out 
