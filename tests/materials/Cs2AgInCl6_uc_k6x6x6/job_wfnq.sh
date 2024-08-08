#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 pw.x     < wfnq.in &> wfnq.in.out 
