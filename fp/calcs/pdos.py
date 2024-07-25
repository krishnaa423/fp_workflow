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
class Pdos:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_pdos = \
f'''&PROJWFC
outdir='./tmp'
prefix='struct'
filpdos='struct_pdos.dat'
/
'''
        
        self.job_pdos = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dos.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dos.job_desc)}projwfc.x -pd .true. < pdos.in &> pdos.in.out 
'''
        
        self.jobs = [
            'job_pdos.sh',
        ]

    def create(self):
        write_str_2_f('pdos.in', self.input_pdos)
        write_str_2_f('job_pdos.sh', self.job_pdos)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_pdos.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'pdos.in*',
            'job_pdos.sh',
            'struct_pdos.dat*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf pdos.in')
        os.system('rm -rf job_pdos.sh')
        
        os.system('rm -rf struct_pdos.dat*')
        os.system('rm -rf pdos.in.out')

#endregion
