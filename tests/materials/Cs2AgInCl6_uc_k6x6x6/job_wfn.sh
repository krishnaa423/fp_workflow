#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 pw.x     < wfn.in &> wfn.in.out 
