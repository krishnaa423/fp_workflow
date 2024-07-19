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
class Kpdos:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_kpdos = \
f'''&PROJWFC
outdir='./tmp'
prefix='struct'
kresolveddos=.true.
filpdos='struct_kpdos.dat'
/
'''
        
        self.job_kpdos = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.kpdos.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.kpdos.job_desc)}projwfc.x -pd .true. < kpdos.in &> kpdos.in.out 
'''
        
        self.jobs = [
            'job_kpdos.sh',
        ]

    def create(self):
        write_str_2_f('kpdos.in', self.input_kpdos)
        write_str_2_f('job_kpdos.sh', self.job_kpdos)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_kpdos.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'kpdos.in*',
            'job_kpdos.sh',
            'struct_kpdos.dat*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf kpdos.in')
        os.system('rm -rf job_kpdos.sh')
        
        os.system('rm -rf struct_kpdos.dat*')
        os.system('rm -rf kpdos.in.out')
#endregion
