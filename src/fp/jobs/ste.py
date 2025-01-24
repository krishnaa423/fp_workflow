#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class SteJob:
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
        if isinstance(self.input_dict['ste']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['ste']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['ste']['job_info'])
    
    def set_inputs_str(self):
        self.input_ste: str = f'''
#region modules
from xctpol.ste import Ste
#endregion

#region variables
#endregion

#region functions
def main():
    ste: Ste = Ste(temp={self.input_dict['ste']['temp']})
    ste.run()
    ste.write()
#endregion

#region classes
#endregion

#region main
main()
#endregion
'''

    def set_jobs_str(self):
        self.job_ste = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

rm -rf ste.out
touch ste.out
exec &> ste.out

{self.scheduler.get_sched_mpi_infix(self.job_info)}python3 script_ste.py &> script_ste.py.out
'''

        self.jobs = [
            './job_ste.sh',
        ]

    def create(self):
        write_str_2_f('script_ste.py', self.input_ste)
        write_str_2_f('job_ste.sh', self.job_ste)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'job_ste.sh',
            'ste.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'script_ste.py',
            'job_ste.sh',

            'ste.h5',
            'ste.out',
            'script_ste.py.out'
        ] 

        for inode in inodes:
            os.system(f'rm -rf ./{inode}')

#endregion
