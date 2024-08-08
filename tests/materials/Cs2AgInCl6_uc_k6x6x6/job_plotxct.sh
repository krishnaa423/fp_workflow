#!/bin/bash




ln -sf WFN_parabands.h5 WFN_fi.h5 
ln -sf WFN_parabands.h5 WFNq_fi.h5 
srun --ntasks=16   --gpus-per-task=1 plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
