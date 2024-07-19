#region: Modules.
from fp.inputs import *
from fp.io import *
from fp.flows import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Inteqp:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_inteqp = \
f'''number_val_bands_coarse {int(self.input.inteqp.val_bands_coarse)}
number_cond_bands_coarse {int(self.input.inteqp.cond_bands_coarse)}
degeneracy_check_override

number_val_bands_fine {int(self.input.inteqp.val_bands_fine)}
number_cond_bands_fine {int(self.input.inteqp.cond_bands_fine)}

use_symmetries_coarse_grid
no_symmetries_fine_grid
'''
        
        self.job_inteqp = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.inteqp.job_desc)}

ln -sf {self.input.inteqp.wfn_co_link} ./WFN_co 
ln -sf {self.input.inteqp.wfn_fi_link} ./WFN_fi 
ln -sf ./eqp1.dat ./eqp_co.dat 
{self.input.scheduler.get_sched_mpi_prefix(self.input.inteqp.job_desc)}inteqp.cplx.x &> inteqp.inp.out 
mv bandstructure.dat bandstructure_inteqp.dat 
'''

        self.jobs = [
            'job_inteqp.sh',   
        ]

    def create(self):
        write_str_2_f('inteqp.inp', self.input_inteqp)
        write_str_2_f('job_inteqp.sh', self.job_inteqp)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_inteqp.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'inteqp.inp*',
            'bandstructure_inteqp.dat',
            'job_inteqp.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf inteqp.inp')
        os.system('rm -rf inteqp.inp.out')
        
        os.system('rm -rf WFN_co')
        os.system('rm -rf WFN_fi')
        os.system('rm -rf eqp_co.dat')
        
        os.system('rm -rf bandstructure_inteqp.dat')
        os.system('rm -rf eqp.dat')
        os.system('rm -rf eqp_q.dat')
        os.system('rm -rf dvmat_norm.dat')
        os.system('rm -rf dcmat_norm.dat')
        
        os.system('rm -rf job_inteqp.sh')
#endregion
