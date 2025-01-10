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
class AbacusJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.abacus_dict_list: List[dict] = [] 
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()
    
    def set_job_info(self):
        if isinstance(self.input_dict['abacus']['common']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['abacus']['common']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['abacus']['common']['job_info'])

    def set_additional_args(self, abacus_dict: dict, calc_dict: dict):
        if calc_dict['args_type']=='extra':
            #Input.
            if calc_dict.get('args', {}).get('input') is not None:
                abacus_dict['input'].update(calc_dict['args']['input'])
            #Input.
            if calc_dict.get('args', {}).get('input') is not None:
                abacus_dict['input'].update(calc_dict['args']['input'])
            #Input.
            if calc_dict.get('args', {}).get('input') is not None:
                abacus_dict['input'].update(calc_dict['args']['input'])
        elif calc_dict['args_type']=='override':
            abacus_dict = calc_dict['args']

        return abacus_dict

    def set_inputs_str(self):
        self.abacus_input_list: List[str] = []
        self.abacus_stru: str = None
        self.abacus_kpt: str = None 
        for calc_idx, calc_dict in enumerate(self.input_dict['abacus']['calculations']):
            # Base.
            abacus_dict  = {
                'input': {
                    'suffix': 'struct',
                    'ntype': self.input.atoms.get_ntyp(),
                    'pseudo_dir': './abacus_pseudos',
                    'orbital_dir': './abacus_orbitals',
                    'ecutwfc': self.input_dict['abacus']['common']['cutoff'],
                    'scf_thr': self.input_dict['abacus']['common']['scf_threshold'],
                    'basis_type': calc_dict['basis'],
                    'calculation': 'scf',
                },
                'structure': {
                    'atomic_species': self.input.atoms.get_abacus_atomic_species(),
                    'numerical_orbital': self.input.atoms.get_abacus_orbitals(),
                    'lattice_constant': [[1.8897259886]],       # value of 1 angstrom in bohr.
                    'lattice_vectors': self.input.atoms.get_abacus_cell(),
                    'atomic_positions': self.input.atoms.get_abacus_atomic_positions(),
                },
                'kpts': self.input_dict['abacus']['common']['kdim']
            }
            
            # Update.
            if self.input_dict['scf']['is_spinorbit']:
                abacus_dict['input']['noncolin'] = 'True'
                abacus_dict['input']['lspinorb'] = 'True'
            abacus_dict = self.set_additional_args(abacus_dict, calc_dict)

            # Write.
            abacus_generator: AbacusInputFile = AbacusInputFile(abacus_dict)
            self.abacus_input_list.append(abacus_generator.get_input_str())
            self.abacus_stru = abacus_generator.get_stru_str()
            self.abacus_kpt = abacus_generator.get_kpt_str()
            abacus_dict = self.abacus_dict_list.append(abacus_dict)

    def set_jobs_str(self):
        abacus_run_string = ''

        for input_idx in range(len(self.abacus_input_list)):
            abacus_run_string += f'cp INPUT{input_idx} INPUT\n'
            abacus_run_string += f'{self.scheduler.get_sched_mpi_prefix(self.job_info)}abacus &> abacus.out\n'

        self.job_abacus: str = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}


{abacus_run_string}
'''

        self.jobs = [
            './job_abacus.sh',
        ]

    def create(self):
        AbacusInputFile(self.abacus_dict_list[0]).copy_pseudos_and_orbitals()
        write_str_2_f('STRU', self.abacus_stru)
        write_str_2_f('KPT', self.abacus_kpt)
        for input_idx, abacus_input in enumerate(self.abacus_input_list):
            write_str_2_f(f'INPUT{input_idx}', abacus_input)
        write_str_2_f('./job_abacus.sh', self.job_abacus)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'INPUT',
            'STRU',
            'KPT',
            'job_abacus.sh',
            'OUT.struct',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf abacus_pseudos')
        os.system('rm -rf abacus_orbitals')
        os.system('rm -rf INPUT*')
        os.system('rm -rf STRU')
        os.system('rm -rf KPT')
        os.system('rm -rf job_abacus.sh')

        os.system('rm -rf abacus.out')
        os.system('rm -rf OUT.*')
        os.system('rm -rf time.json')

#endregion
