#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
from fp.inputs.bgw import BgwInputFile
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPhJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['xctph']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['xctph']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['xctph']['job_info'])
    
    def set_inputs_str(self):
        self.input_xctph: str = '''
#region modules
from xctph.xctph import Xctph
#endregion

#region variables
#endregion

#region functions
def main():
    xctph: Xctph = Xctph()
    xctph.calc()
    xctph.write()
#endregion

#region classes
#endregion

#region main
main()
#endregion
'''

    def set_jobs_str(self):
        self.job_xctph = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

rm -rf xctph.out
touch xctph.out
exec &> xctph.out

{self.scheduler.get_sched_mpi_prefix(self.job_info)}python3 script_xctph.py &> script_xctph.py.out
'''

        self.jobs = [
            './job_xctph.sh',
        ]

    def create(self):
        write_str_2_f('script_xctph.py', self.input_xctph)
        write_str_2_f('job_xctph.sh', self.job_xctph)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'job_xctph.sh',
            'xct.h5',
            'eph*.h5',
            'xctph*.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'script_xctph.py',
            'job_xctph.sh',

            'xct.h5',
            'eph*.h5',
            'xctph*.h5',
            'xctph.out',
            'script_xctph.py.out'
        ] 

        for inode in inodes:
            os.system(f'rm -rf ./{inode}')

#endregion
