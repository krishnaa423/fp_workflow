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
class Sigma:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_sigma = \
f'''# kpoints
{self.input.sigma.get_kgrid_str(self.input.wfn)}
no_symmetries_q_grid
# use_symmetries_q_grid

# Bands.
number_bands {self.input.sigma.bands}
band_index_min {self.input.sigma.band_min}
band_index_max {self.input.sigma.band_max}
degeneracy_check_override

# G-cutoff
screened_coulomb_cutoff {self.input.sigma.cutoff}


# Options
dont_use_vxcdat

# IO
# verbosity 3
use_wfn_hdf5
'''
        
        self.job_sigma = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.sigma.job_desc)}

ln -sf {self.input.sigma.wfn_inner_link} ./WFN_inner.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.sigma.job_desc)}sigma.cplx.x &> sigma.inp.out
'''
        
        self.jobs = [
            'job_sigma.sh',
        ]

    def create(self):
        write_str_2_f('sigma.inp', self.input_sigma)
        write_str_2_f('job_sigma.sh', self.job_sigma)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_sigma.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'sigma.inp*',
            'job_sigma.sh',
            'eqp0.dat',
            'eqp1.dat',
            'sigma_hp.log',
            'ch_converge.dat',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf sigma.inp')
        os.system('rm -rf job_sigma.sh')
        
        os.system('rm -rf ./WFN_inner.h5')
        os.system('rm -rf eqp0.dat')
        os.system('rm -rf eqp1.dat')
        os.system('rm -rf sigma_hp.log')
        os.system('rm -rf ch_converge.dat')
        os.system('rm -rf sigma.inp.out')
        os.system('rm -rf dtmat')
        os.system('rm -rf vxc.dat')
        os.system('rm -rf x.dat')
#endregion
