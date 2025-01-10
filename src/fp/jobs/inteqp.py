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
class InteqpJob:
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
        if isinstance(self.input_dict['inteqp']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['inteqp']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['inteqp']['job_info'])

    def set_inputs_str(self):
        # Base.
        input_inteqp_dict: dict = {
            'maps': {
                'number_val_bands_coarse': self.input_dict['inteqp']['num_val_bands'],
                'number_cond_bands_coarse': self.input_dict['dftelbands']['num_cond_bands'],
                'number_val_bands_fine': self.input_dict['inteqp']['num_val_bands'] - 1,
                'number_cond_bands_fine': self.input_dict['dftelbands']['num_cond_bands'],
                'degeneracy_check_override': '',
                'use_symmetries_coarse_grid': '',
                'no_symmetries_fine_grid': '',
            },
        }

        # Additions.
        args_dict = self.input_dict['inteqp']['args']
        args_type = self.input_dict['inteqp']['args_type']
        input_inteqp_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_inteqp_dict
        )

        # Write.
        self.input_inteqp: str = BgwInputFile.write_general(input_inteqp_dict)

    def set_jobs_str(self):
        self.job_inteqp = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf {self.input_dict['inteqp']['wfnco_link']} ./WFN_co 
ln -sf {self.input_dict['inteqp']['wfnfi_link']} ./WFN_fi 
ln -sf ./eqp1.dat ./eqp_co.dat 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}inteqp.cplx.x &> inteqp.inp.out 
mv bandstructure.dat bandstructure_inteqp.dat 
'''

        self.jobs = [
            './job_inteqp.sh',   
        ]

    def create(self):
        write_str_2_f('inteqp.inp', self.input_inteqp)
        write_str_2_f('job_inteqp.sh', self.job_inteqp)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'inteqp.inp*',
            'bandstructure_inteqp.dat',
            'job_inteqp.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf inteqp.inp')
        os.system('rm -rf inteqp.inp.out')
        
        os.system('rm -rf WFN_co')
        os.system('rm -rf WFN_fi')
        os.system('rm -rf eqp_co.dat')
        
        os.system('rm -rf bandstructure_inteqp.dat')
        os.system('rm -rf eqp.dat')
        os.system('rm -rf eqp_q.dat')
        os.system('rm -rf dvmat_norm.dat')
        os.system('rm -rf dcmat_norm.dat')
        
        os.system('rm -rf job_inteqp.sh')
#endregion
