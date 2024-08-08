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
class Esf:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.script_esf = \
f'''
from fp.analysis.esf import Esf

esf = Esf()
esf.calc()
esf.write()
'''
        
        self.job_esf = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.esf.job_desc)}

python3 script_esf.py &> script_esf.out
'''

        self.jobs = [
            'job_esf.sh',
        ]

    def create(self):
        write_str_2_f('script_esf.py', self.script_esf)
        write_str_2_f('job_esf.sh', self.job_esf)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_esf.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'script_esf.py',
            'job_esf.sh',

            'script_esf.out',
            'esf.h5',
            'esf.xsf',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'script_esf.py',
            'job_esf.sh',

            'script_esf.out',
            'esf.h5',
            'esf.xsf',
        ] 

        for inode in inodes:
            os.system(f'rm -rf ./{inode}')

#endregion
