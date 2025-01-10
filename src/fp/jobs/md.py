#region modules
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
from fp.inputs.abacus import AbacusInputFile
from importlib.util import find_spec
from typing import List 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class MdJob:
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
        if isinstance(self.input_dict['md']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['md']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['md']['job_info'])

    def set_inputs_str(self):
        #Base. 
        input_md_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'md',
                    'tprnfor': True,
                    'dt': self.input_dict['md']['time_step'],
                    'nstep': self.input_dict['md']['number_steps'],
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                    'nosym': True,
                },
                'electrons': {},
                'ions': {
                    'pot_extrapolation': 'second-order',
                    'wfc_extrapolation': 'second-order',
                    'ion_temperature': 'initial',
                    'tempw': self.input_dict['md']['temperature'],
                },
                'cell': {},
            },
            'blocks': {
                'atomic_species': self.input.atoms.get_qe_scf_atomic_species(),
                'cell_parameters': self.input.atoms.get_qe_scf_cell(),
                'atomic_positions': self.input.atoms.get_qe_scf_atomic_positions(),
                'kpoints': [1, 1, 1],           # MD will be gamma point. 
            },
            'kpoints_type': 'automatic',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_md_dict['namelists']['system']['noncolin'] = True
            input_md_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['md']['args']
        args_type = self.input_dict['md']['args_type']
        input_md_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_md_dict
        )

        # Get string. 
        self.input_md: str = QePwInputFile(input_md_dict, self.input_dict).get_input_str()

    def set_jobs_str(self):
        self.job_scf: str = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < md.in &> md.in.out

cp ./tmp/struct.save/data-file-schema.xml ./md.xml
'''
        self.jobs = [
            './job_md.sh',
        ]
        
    def create(self):
        write_str_2_f('md.in', self.input_md)
        write_str_2_f('job_md.sh', self.job_scf)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'md.in',
            'job_md.sh',
            './tmp'
            'md.xml',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf md.in')
        os.system('rm -rf job_md.sh')

        os.system('rm -rf md.in.out')
        os.system('rm -rf md.xml')

#endregion
