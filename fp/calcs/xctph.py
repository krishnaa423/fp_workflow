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
class XctPh:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.script_xctph = \
f'''
from fp.analysis.xctph import XctPh

xctph = XctPh(
    './struct_elph_',
    './bseq', 
    './fullgridflow.pkl',
    './input.pkl',
)
xctph.get_xctph()
xctph.write()
'''
        
        self.job_xctph = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.xctph.job_desc)}

python3 script_xctph.py &> script_xctph.out
'''

        self.jobs = [
            'job_xctph.sh',
        ]

    def create(self):
        write_str_2_f('script_xctph.py', self.script_xctph)
        write_str_2_f('job_xctph.sh', self.job_xctph)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_xctph.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'script_xctph.py',
            'job_xctph.sh',

            'script_xctph.out',
            'xctph.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'script_xctph.py',
            'job_xctph.sh',

            'script_xctph.out',
            'xctph.h5',
        ] 

        for inode in inodes:
            os.system(f'rm -rf ./{inode}')

#endregion
