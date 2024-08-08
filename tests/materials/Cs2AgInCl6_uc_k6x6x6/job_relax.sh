#!/bin/bash



srun --ntasks=16   --gpus-per-task=1 pw.x     < relax.in &> relax.in.out

cp ./tmp/struct.save/data-file-schema.xml ./relax.xml
