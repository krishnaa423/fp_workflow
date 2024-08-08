#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 pw.x     < wfnfi.in &> wfnfi.in.out 
