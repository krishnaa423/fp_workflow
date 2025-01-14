#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPolJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.script_xctpol = \
f'''
from fp.analysis.xctpol import XctPolCalc

xctpol = XctPolCalc(
    './eph_xctph.h5',
    './xctph_elhole.h5', 
    './xctph_hole.h5',
    max_error={self.input.xctpol.max_error},
    max_steps={self.input.xctpol.max_steps},
)
xctpol.get_xctpol()
xctpol.get_xctpol_energy()
xctpol.write()
'''
        
        self.job_xctpol = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

python3 script_xctpol.py &> script_xctpol.out
'''

        self.jobs = [
            './job_xctpol.sh',
        ]

    def create(self):
        write_str_2_f('script_xctpol.py', self.script_xctpol)
        write_str_2_f('job_xctpol.sh', self.job_xctpol)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'script_xctpol.py',
            'job_xctpol.sh',

            'script_xctpol.out',
            'xctpol.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'script_xctpol.py',
            'job_xctpol.sh',

            'script_xctpol.out',
            'xctpol.h5',
        ] 

        for inode in inodes:
            os.system(f'rm -rf ./{inode}')

#endregion
