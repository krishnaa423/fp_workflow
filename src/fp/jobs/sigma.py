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
class SigmaJob:
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
        if isinstance(self.input_dict['sig']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['sig']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['sig']['job_info'])

    def set_inputs_str(self):
        # Base.
        input_sigma_dict: dict = {
            'maps': {
                'no_symmetries_q_grid': '',
                'number_bands': self.input_dict['sig']['total_cond_bands'] + self.input_dict['total_valence_bands'],
                'band_index_min': self.input_dict['total_valence_bands'] - self.input_dict['sig']['se_val_bands'] + 1,
                'band_index_max': self.input_dict['total_valence_bands'] + self.input_dict['sig']['se_cond_bands'],
                'degeneracy_check_override': '',
                'screened_coulomb_cutoff': self.input_dict['sig']['cutoff'],
                'dont_use_vxcdat': '',
                'use_wfn_hdf5': '',
            },
            'blocks': {
                'kpoints': self.input.sig.get_kgrid(self.input.wfn)
            }
        }

        # Additions.
        args_dict = self.input_dict['sig']['args']
        args_type = self.input_dict['sig']['args_type']
        input_sigma_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_sigma_dict
        )

        # Write.
        self.input_sigma: str = BgwInputFile.write_general(input_sigma_dict)

    def set_jobs_str(self):
        self.job_sigma = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf {self.input_dict['sig']['wfninner_link']} ./WFN_inner.h5 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}sigma.cplx.x &> sigma.inp.out
'''
        
        self.jobs = [
            './job_sigma.sh',
        ]

    def create(self):
        write_str_2_f('sigma.inp', self.input_sigma)
        write_str_2_f('job_sigma.sh', self.job_sigma)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'sigma.inp*',
            'job_sigma.sh',
            'eqp0.dat',
            'eqp1.dat',
            'sigma_hp.log',
            'ch_converge.dat',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf sigma.inp')
        os.system('rm -rf job_sigma.sh')
        
        os.system('rm -rf ./WFN_inner.h5')
        os.system('rm -rf eqp0.dat')
        os.system('rm -rf eqp1.dat')
        os.system('rm -rf sigma_hp.log')
        os.system('rm -rf ch_converge.dat')
        os.system('rm -rf sigma.inp.out')
        os.system('rm -rf dtmat')
        os.system('rm -rf vxc.dat')
        os.system('rm -rf x.dat')
#endregion
