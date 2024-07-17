#!/bin/bash


mpirun -n 1 pw.x     < relax.in &> relax.in.out

cp ./tmp/struct.save/data-file-schema.xml ./relax.xml
