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
class KernelJob:
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
        if isinstance(self.input_dict['ker']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['ker']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['ker']['job_info'])

    def set_inputs_str(self):
        # Base.
        input_kernel_dict: dict = {
            'maps': {
                'exciton_Q_shift': [
                    2,
                    self.input_dict['abs']['Qshift'][0],
                    self.input_dict['abs']['Qshift'][1],
                    self.input_dict['abs']['Qshift'][2],
                ],
                'use_symmetries_coarse_grid': '',
                'number_val_bands': self.input_dict['abs']['num_val_bands'],
                'number_cond_bands': self.input_dict['abs']['num_cond_bands'],
                'use_wfn_hdf5': '',
            },
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_kernel_dict['maps']['spinor'] = ''
        args_dict = self.input_dict['ker']['args']
        args_type = self.input_dict['ker']['args_type']
        input_kernel_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_kernel_dict
        )

        # Write.
        self.input_kernel: str = BgwInputFile.write_general(input_kernel_dict)

    def set_jobs_str(self):
        self.job_kernel = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf {self.input_dict['abs']['wfnco_link']} WFN_co.h5
ln -sf {self.input_dict['abs']['wfnqco_link']} WFNq_co.h5
{self.scheduler.get_sched_mpi_prefix(self.job_info)}kernel.cplx.x &> kernel.inp.out
'''

        self.jobs = [
            './job_kernel.sh',
        ]

    def create(self):
        write_str_2_f('kernel.inp', self.input_kernel)
        write_str_2_f('job_kernel.sh', self.job_kernel)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'kernel.inp*',
            'bsemat.h5',
            'job_kernel.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf kernel.inp')
        os.system('rm -rf job_kernel.sh')
        
        os.system('rm -rf ./WFN_co.h5')
        os.system('rm -rf bsemat.h5')
        os.system('rm -rf kernel.inp.out')
#endregion
