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
class PhmodesJob:
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
        if isinstance(self.input_dict['phmodes']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['phmodes']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['phmodes']['job_info'])
    
    def set_inputs_str(self):
        # Base.
        phmodes_qidx = self.input_dict['phmodes']['qpt_idx']
        input_dynmat_dict = {
            'namelists': {
                'input': {
                    'asr': "'crystal'",
                    'fildyn': f"'struct.dyn{phmodes_qidx}'",
                    'filxsf': "'struct_phmodes.axsf'",
                }
            }
        }

        # Additions.
        args_dict = self.input_dict['phmodes']['args']
        args_type = self.input_dict['phmodes']['args_type']
        input_dynmat_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_dynmat_dict
        )

        # Write.
        self.input_dynmat: str = QePwInputFile.write_general(input_dynmat_dict)

    def set_jobs_str(self):
        self.job_dynmat = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}dynmat.x < dynmat.in &> dynmat.in.out 
'''
        self.jobs = [
            './job_dynmat.sh',
        ]

    def create(self):
        write_str_2_f('dynmat.in', self.input_dynmat)
        write_str_2_f('job_dynmat.sh', self.job_dynmat)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'dynmat.in*',
            'job_dynmat.sh',
            'struct.dyn*',
            'struct_phmodes.axsf',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf dynmat.in')
        os.system('rm -rf dynmat.out')
        os.system('rm -rf dynmat.mold')
        os.system('rm -rf input_tmp.in')
        os.system('rm -rf dynmat.in.out')
        os.system('rm -rf job_dynmat.sh')
        
        os.system('rm -rf struct.dyn*')
        os.system('rm -rf struct_phmodes.axsf')

#endregion
