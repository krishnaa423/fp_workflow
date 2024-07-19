#region: Modules.
from fp.inputs import *
from fp.io import *
from fp.flows import *
from pkg_resources import resource_filename
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Phmodes:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_dynmat = \
f'''&INPUT
asr='crystal'
fildyn='struct.dyn{self.input.phmodes.qidx}'
filxsf='struct_phmodes.axsf'
/
'''
        
        self.job_dynmat = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.phmodes.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.phmodes.job_desc)}dynmat.x < dynmat.in &> dynmat.in.out 
'''
        self.jobs = [
            'job_dynmat.sh',
        ]

    def create(self):
        write_str_2_f('dynmat.in', self.input_dynmat)
        write_str_2_f('job_dynmat.sh', self.job_dynmat)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_dynmat.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'dynmat.in*',
            'job_dynmat.sh',
            'struct.dyn*',
            'struct_phmodes.axsf',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf dynmat.in')
        os.system('rm -rf dynmat.out')
        os.system('rm -rf dynmat.mold')
        os.system('rm -rf input_tmp.in')
        os.system('rm -rf dynmat.in.out')
        os.system('rm -rf job_dynmat.sh')
        
        os.system('rm -rf struct.dyn*')
        os.system('rm -rf struct_phmodes.axsf')
#endregion
