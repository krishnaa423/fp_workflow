#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
from fp.inputs.wannier import WannierWinFile
from importlib.util import find_spec
from typing import List 
from fp.structure.kpts import Kgrid
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class WannierJob:
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
        if isinstance(self.input_dict['wannier']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['wannier']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['wannier']['job_info'])

    def update_args(self):
        # wfnwan.
        if self.input_dict.get('wannier', {}).get('args', {}).get('wfnwan') is not None:
            match self.input_dict['wannier']['args_type']:
                case 'extra':
                    if self.input_dict['wannier']['args']['wfnwan'].get('namelists') is not None:
                        for key, value in self.input_dict['wannier']['args']['wfnwan']['namelists'].items():
                            self.input_wfnwan_dict['namelists'][key].update(value)
                    if self.input_dict['wannier']['args']['wfnwan'].get('blocks') is not None:
                        self.input_wfnwan_dict['blocks'].update(self.input_dict['wannier']['args']['wfnwan']['blocks'])
                case 'override':
                    self.input_wfnwan_dict = self.input_dict['wannier']['args']['wfnwan']

        # wan.
        if self.input_dict.get('wannier', {}).get('args', {}).get('wan') is not None:
            match self.input_dict['wannier']['args_type']:
                case 'extra':
                    if self.input_dict['wannier']['args']['wan'].get('maps') is not None:
                        self.input_wan_dict['maps'].update(self.input_dict['wannier']['args']['wan']['maps'])
                    if self.input_dict['wannier']['args']['wan'].get('blocks') is not None:
                        self.input_wan_dict['blocks'].update(self.input_dict['wannier']['args']['wan']['blocks'])
                case 'override':
                    self.input_wan_dict = self.input_dict['wannier']['args']['wan']

        # pw2wan.
        if self.input_dict.get('wannier', {}).get('args', {}).get('pw2wan') is not None:
            match self.input_dict['wannier']['args_type']:
                case 'extra':
                    if self.input_dict['wannier']['args']['pw2wan'].get('namelists') is not None:
                        for key, value in self.input_dict['wannier']['args']['pw2wan']['namelists'].items():
                            self.input_pw2wan_dict['namelists'][key].update(value)
                case 'override':
                    self.input_pw2wan_dict = self.input_dict['wannier']['args']['pw2wan']

    def get_exclude_bands(self):
        exclude_bands = None
        total_valence_bands = self.input_dict['total_valence_bands']
        num_cond_bands = self.input_dict['wannier']['num_cond_bands']
        num_val_bands = self.input_dict['wannier']['num_val_bands']

        # We will only exclude valence bands, since conduction bands in wfnwan.in and wan.win are set to  
        #always match. 
        if total_valence_bands > num_val_bands:
            exclude_bands = [
                (1, total_valence_bands - num_val_bands),
            ]

        return exclude_bands

    def set_inputs_str(self):
        #Base. 
        kpts = Kgrid(self.input.atoms, self.input_dict['wannier']['kdim']).get_kpts()
        self.input_wfnwan_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'bands',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'nbnd': self.input_dict['total_valence_bands'] + self.input_dict['wannier']['num_cond_bands'],
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                },
                'electrons': {},
                'ions': {},
                'cell': {},
            },
            'blocks': {
                'atomic_species': self.input.atoms.get_qe_scf_atomic_species(),
                'cell_parameters': self.input.atoms.get_qe_scf_cell(),
                'atomic_positions': self.input.atoms.get_qe_scf_atomic_positions(),
                'kpoints': kpts,
            },
            'kpoints_type': 'crystal',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        self.input_wan_dict: dict = {
            'maps': {
                'mp_grid': self.input_dict['wannier']['kdim'],
                'num_bands': self.input_dict['wannier']['num_val_bands']+ self.input_dict['wannier']['num_cond_bands'],
                'exclude_bands': self.get_exclude_bands(),
                'num_wann': self.input_dict['wannier']['num_val_bands'] + self.input_dict['wannier']['num_cond_bands'],
                'auto_projections': '.true.',
                'wannier_plot': '.true.',
                'write_hr': '.true.',
                'write_u_matrices': '.true.',
            },
            'blocks': {
                'unit_cell_cart': self.input.atoms.get_wan_cell(),
                'atoms_cart': self.input.atoms.get_wan_atomic_positions(),
                'kpoints': kpts,
            }
        }
        self.input_pw2wan_dict: dict = {
            'namelists': {
                'inputpp': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'seedname': "'wan'",
                    'write_amn': '.true.',
                    'write_mmn': '.true.',
                    'write_unk': '.true.',
                    'scdm_proj': '.true.',
                }
            }
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            self.input_wfnwan_dict['namelists']['system']['noncolin'] = True
            self.input_wfnwan_dict['namelists']['system']['lspinorb'] = True
            self.input_wan_dict['maps']['spinors'] = '.true.'
        #override or extra. 
        self.update_args()

        # Get string. 
        self.input_wfnwan: str = QePwInputFile(self.input_wfnwan_dict, self.input_dict).get_input_str()
        self.input_wan: str = WannierWinFile.write_general(self.input_wan_dict)
        self.input_pw2wan: str = QePwInputFile.write_general(self.input_pw2wan_dict)

    def set_jobs_str(self):
        self.job_wfnwan = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x < wfnwan.in &> wfnwan.in.out 
'''
        
        self.job_wanpp = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}wannier90.x {self.scheduler.get_sched_mpi_infix(self.job_info)} -pp wan &> wan.win.pp.out
'''
        
        self.job_pw2wan = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw2wannier90.x < pw2wan.in &> pw2wan.in.out 
'''
        
        self.job_wan = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}wannier90.x {self.scheduler.get_sched_mpi_infix(self.job_info)} wan  &> wan.win.out 
'''
        
        self.jobs = [
            './job_wfnwan.sh',
            './job_wanpp.sh',
            './job_pw2wan.sh',
            './job_wan.sh',
        ]

    def create(self):
        write_str_2_f('wfnwan.in', self.input_wfnwan)
        write_str_2_f('job_wfnwan.sh', self.job_wfnwan)
        write_str_2_f('wan.win', self.input_wan)
        write_str_2_f('job_wanpp.sh', self.job_wanpp)
        write_str_2_f('pw2wan.in', self.input_pw2wan)
        write_str_2_f('job_pw2wan.sh', self.job_pw2wan)
        write_str_2_f('job_wan.sh', self.job_wan)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wan*',
            'wfnwan.in*',
            'job_wfnwan.sh',
            'job_wanpp.sh',
            'job_pw2wan.sh',
            'job_wan.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wan*')
        os.system('rm -rf UNK*')
        
        os.system('rm -rf pw2wan.in')
        os.system('rm -rf pw2wan.in.out')
        os.system('rm -rf wfnwan.in')
        os.system('rm -rf wfnwan.in.out')
        
        
        os.system('rm -rf job_wfnwan.sh')
        os.system('rm -rf job_wanpp.sh')
        os.system('rm -rf job_pw2wan.sh')
        os.system('rm -rf job_wan.sh')

#endregion
