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
class EpsilonJob:
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
        if isinstance(self.input_dict['eps']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['eps']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['eps']['job_info'])

    def set_inputs_str(self):
        # Base.
        input_epsilon_dict: dict = {
            'maps': {
                'number_bands': self.input_dict['eps']['num_cond_bands'] + self.input_dict['total_valence_bands'],
                'degeneracy_check_override': '',
                'epsilon_cutoff': self.input_dict['eps']['cutoff'],
                'use_wfn_hdf5': '',
            },
            'blocks': {
                'qpoints': self.input.eps.get_qgrid_str(self.input.wfn, self.input_dict['wfnq']['qshift'])
            }
        }

        # Additions.
        args_dict = self.input_dict['eps']['args']
        args_type = self.input_dict['eps']['args_type']
        input_epsilon_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_epsilon_dict
        )

        # Write.
        self.input_epsilon: str = BgwInputFile.write_general(input_epsilon_dict)

    def set_jobs_str(self):
        self.job_epsilon = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf {self.input_dict['eps']['wfnlink']} ./WFN.h5 
ln -sf {self.input_dict['eps']['wfnqlink']} ./WFNq.h5 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}epsilon.cplx.x &> epsilon.inp.out 
'''

        self.jobs = [
            './job_epsilon.sh',
        ]

    def create(self):
        write_str_2_f('epsilon.inp', self.input_epsilon)
        write_str_2_f('job_epsilon.sh', self.job_epsilon)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'epsilon.inp*',
            'job_epsilon.sh',
            'epsmat.h5',
            'eps0mat.h5',
            'epsilon.log',
            'chi_converge.dat',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf epsilon.inp')
        os.system('rm -rf job_epsilon.sh')
        
        os.system('rm -rf ./WFN.h5')
        os.system('rm -rf ./WFNq.h5')
        os.system('rm -rf ./epsmat.h5')
        os.system('rm -rf ./eps0mat.h5')
        os.system('rm -rf epsilon.log')
        os.system('rm -rf chi_converge.dat')
        os.system('rm -rf epsilon.inp.out')
        
        os.system('rm -rf checkbz.log')

#endregion
