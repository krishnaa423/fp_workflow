#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 pw.x     < wfnqfi.in &> wfnqfi.in.out 
