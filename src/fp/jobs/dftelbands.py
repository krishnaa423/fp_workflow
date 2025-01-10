#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
from fp.structure.kpath import KPath
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DftelbandsJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.job_pw2bgw_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['dftelbands']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['dftelbands']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['dftelbands']['job_info'])

        # pw2bgw.
        if isinstance(self.input_dict['dftelbands']['job_pw2bgw_info'], str):
            self.job_pw2bgw_info = JobProcDesc.from_job_id(
                self.input_dict['dftelbands']['job_pw2bgw_info'],
                self.input_dict,
            )
        else:
            self.job_pw2bgw_info = JobProcDesc(**self.input_dict['dftelbands']['job_pw2bgw_info'])

    def set_inputs_str(self):
        input_dftelbands_dict = {
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
                    'nbnd': self.input_dict['dftelbands']['num_cond_bands'] + self.input_dict['total_valence_bands'],
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
                'kpoints': self.input.dftelbands.get_kgrid(kpath=KPath(self.input_dict, self.input.atoms.atoms)),
            },
            'kpoints_type': 'crystal_b',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }

        input_pw2bgw_dict = {
            'namelists': {
                'input_pw2bgw': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'real_or_complex': '2',
                    'wfng_flag': '.true.',
                    'wfng_file': "'WFN_dftelbands'",
                    'wfng_kgrid': '.true.',
                    'wfng_nk1': 0,
                    'wfng_nk2': 0,
                    'wfng_nk3': 0,
                    'wfng_dk1': 0.0,
                    'wfng_dk2': 0.0,
                    'wfng_dk3': 0.0,
                }
            }
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_dftelbands_dict['namelists']['system']['noncolin'] = True
            input_dftelbands_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['dftelbands']['args']
        args_type = self.input_dict['dftelbands']['args_type']
        input_dftelbands_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_dftelbands_dict
        )
        args_dict = self.input_dict['dftelbands']['pw2bgw_args']
        args_type = self.input_dict['dftelbands']['pw2bgw_args_type']
        input_pw2bgw_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_pw2bgw_dict
        )

        # Get string. 
        self.input_dftelbands: str = QePwInputFile(input_dftelbands_dict, self.input_dict).get_input_str()
        self.input_pw2bgw: str = QePwInputFile.write_general(input_pw2bgw_dict)

    def set_jobs_str(self):
        self.job_dftelbands = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x < dftelbands.in &> dftelbands.in.out  
cp ./tmp/struct.xml ./dftelbands.xml 
'''     
        self.job_dftelbands_pw2bgw = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_pw2bgw_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_pw2bgw_info)}pw2bgw.x -pd .true. < dftelbands_pw2bgw.in &> dftelbands_pw2bgw.in.out 
cp ./tmp/WFN_dftelbands ./
wfn2hdf.x BIN WFN_dftelbands WFN_dftelbands.h5
'''
        
        self.jobs = [
            './job_dftelbands.sh',
            './job_dftelbands_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f('dftelbands.in', self.input_dftelbands)
        write_str_2_f('job_dftelbands.sh', self.job_dftelbands)
        write_str_2_f('dftelbands_pw2bgw.in', self.input_pw2bgw)
        write_str_2_f('job_dftelbands_pw2bgw.sh', self.job_dftelbands_pw2bgw)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'dftelbands.in*',
            'job_dftelbands.sh',
            'dftelbands_pw2bgw.in*',
            'job_dftelbands_pw2bgw.sh',
            'WFN_dftelbands',
            'WFN_dftelbands.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf dftelbands.in')
        os.system('rm -rf job_dftelbands.sh')
        os.system('rm -rf dftelbands_pw2bgw.in')
        os.system('rm -rf job_dftelbands_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf dftelbands.in.out')
        os.system('rm -rf dftelbands_pw2bgw.in.out')
        os.system('rm -rf dftelbands.xml')
        os.system('rm -rf kgrid.inp kgrid.log kgrid.out')
        os.system('rm -rf WFN_dftelbands')
        os.system('rm -rf WFN_dftelbands.h5')
#endregion
