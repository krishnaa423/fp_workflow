#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PdosJob:
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
        if isinstance(self.input_dict['pdos']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['pdos']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['pdos']['job_info'])

    def set_inputs_str(self):
        #Base.
        input_pdos_dict = {
            'namelists': {
                'projwfc': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'filpdos': "'struct_pdos.dat'",
                }
            }
        }

        #Additions.
        #override or extra. 
        args_dict = self.input_dict['pdos']['args']
        args_type = self.input_dict['pdos']['args_type']
        input_pdos_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_pdos_dict
        )

        #Write.
        self.input_pdos: str = QePwInputFile.write_general(input_pdos_dict)

    def set_jobs_str(self):
        self.job_pdos = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}projwfc.x -pd .true. < pdos.in &> pdos.in.out 
'''
        
        self.jobs = [
            './job_pdos.sh',
        ]

    def create(self):
        write_str_2_f('pdos.in', self.input_pdos)
        write_str_2_f('job_pdos.sh', self.job_pdos)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
