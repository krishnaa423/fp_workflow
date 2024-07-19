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
class Epsilon:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_epsilon = \
f'''# Qpoints 
{self.input.epsilon.get_qgrid_str(self.input.wfn, self.input.wfnq.qshift)}

# Bands
number_bands {self.input.epsilon.bands}
degeneracy_check_override

# G-Cutoff. 
epsilon_cutoff {self.input.epsilon.cutoff}

# Options

# IO. 
use_wfn_hdf5
'''
        
        self.job_epsilon = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.epsilon.job_desc)}

ln -sf {self.input.epsilon.wfn_link} ./WFN.h5 
ln -sf {self.input.epsilon.wfnq_link} ./WFNq.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.epsilon.job_desc)}epsilon.cplx.x &> epsilon.inp.out 
'''

        self.jobs = [
            'job_epsilon.sh',
        ]

    def create(self):
        write_str_2_f('epsilon.inp', self.input_epsilon)
        write_str_2_f('job_epsilon.sh', self.job_epsilon)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_epsilon.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'epsilon.inp*',
            'job_epsilon.sh',
            'epsmat.h5',
            'eps0mat.h5',
            'epsilon.log',
            'chi_converge.dat',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf epsilon.inp')
        os.system('rm -rf job_epsilon.sh')
        
        os.system('rm -rf ./WFN.h5')
        os.system('rm -rf ./WFNq.h5')
        os.system('rm -rf ./epsmat.h5')
        os.system('rm -rf ./eps0mat.h5')
        os.system('rm -rf epsilon.log')
        os.system('rm -rf chi_converge.dat')
        os.system('rm -rf epsilon.inp.out')
        
        os.system('rm -rf checkbz.log')
#endregion
