#!/bin/bash


mpirun -n 1 pw.x     < scf.in &> scf.in.out

cp ./tmp/struct.save/data-file-schema.xml ./scf.xml
