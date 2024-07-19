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
class Kernel:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_kernel = \
f'''# Q-points
exciton_Q_shift 2 {self.input.kernel.Qshift[0]:15.10f} {self.input.kernel.Qshift[1]:15.10f} {self.input.kernel.Qshift[2]:15.10f}
use_symmetries_coarse_grid

# Bands 
number_val_bands {self.input.kernel.val_bands_coarse}
number_cond_bands {self.input.kernel.cond_bands_coarse}
#spinor

# Options
#extended_kernel

# IO. 
use_wfn_hdf5
'''
        
        self.job_kernel = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.kernel.job_desc)}

ln -sf {self.input.kernel.wfn_co_link} WFN_co.h5
ln -sf {self.input.kernel.wfnq_co_link} WFNq_co.h5
{self.input.scheduler.get_sched_mpi_prefix(self.input.kernel.job_desc)}kernel.cplx.x &> kernel.inp.out
'''

        self.jobs = [
            'job_kernel.sh',
        ]

    def create(self):
        write_str_2_f('kernel.inp', self.input_kernel)
        write_str_2_f('job_kernel.sh', self.job_kernel)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_kernel.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'kernel.inp*',
            'bsemat.h5',
            'job_kernel.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf kernel.inp')
        os.system('rm -rf job_kernel.sh')
        
        os.system('rm -rf ./WFN_co.h5')
        os.system('rm -rf bsemat.h5')
        os.system('rm -rf kernel.inp.out')
#endregion
