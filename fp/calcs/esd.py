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
class EsdStep:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.script_esd_step = \
f'''
from fp.analysis.esd import Esd
from ase.io import Atoms 

Esd.static_step(fmax={self.input.esdstep.fmax}, max_steps={self.input.esdstep.max_steps})
'''
        
        self.job_esd_step = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.esdstep.job_desc)}

python3 script_esd_step.py &> script_esd_step.out
'''
        
        self.jobs = [
            'job_esd_step.sh',
        ]

    def create(self):
        write_str_2_f('script_esd_step.py', self.script_esd_step)
        write_str_2_f('job_esd_step.sh', self.job_esd_step)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_esd_step.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'script_esd_step*',
            'job_esd_step.sh',
            'esd_atoms.xsf',
            'esd.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'script_esd_step*',
            'job_esd_step.sh',
            'esd_atoms.xsf',
            'esd.h5',
        ] 

        for inode in inodes:
            os.system(f'rm -rf {inode}')

#endregion
