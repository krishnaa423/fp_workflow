#region modules
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class ElphJob:
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
        if isinstance(self.input_dict['elph']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['elph']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['elph']['job_info'])

    def set_inputs_str(self):
        self.input_elph: str = '''
#region modules
from xctph.elph import Elph
#endregion

#region variables
#endregion

#region functions
def main():
    elph: Elph = Elph()
    elph.read()
    elph.write()
#endregion

#region classes
#endregion

#region main
main()
#endregion
'''

    def set_jobs_str(self):
        self.job_elph = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

rm -rf elph.out
touch elph.out
exec &> elph.out

{self.scheduler.get_sched_mpi_prefix(self.job_info)}python3 script_elph.py &> script_elph.py.out 
'''
        
        self.jobs = [
            './job_elph.sh',
        ]

    def create(self):
        write_str_2_f('script_elph.py', self.input_elph)
        write_str_2_f('job_elph.sh', self.job_elph)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'elph.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf ./script_elph.py')
        os.system('rm -rf ./job_elph.sh')

        os.system('rm -rf ./script_elph.py.out')
        os.system('rm -rf ./elph.h5')
        os.system('rm -rf ./elph.out')
#endregion