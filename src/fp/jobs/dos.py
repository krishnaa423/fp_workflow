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
class DosJob:
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
        if isinstance(self.input_dict['dos']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['dos']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['dos']['job_info'])

    def set_inputs_str(self):
        #Base.
        input_wfndos_dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'bands',
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'nbnd': self.input_dict['dos']['num_cond_bands'] + self.input_dict['total_valence_bands'],
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
                'kpoints': self.input_dict['dos']['kdim'],
            },
            'kpoints_type': 'automatic',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        input_dos_dict = {
            'namelists': {
                'dos': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'fildos': "'struct_dos.dat'",
                }
            }
        }

        #Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_wfndos_dict['namelists']['system']['noncolin'] = True
            input_wfndos_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['dos']['wfn_args']
        args_type = self.input_dict['dos']['wfn_args_type']
        input_wfndos_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_wfndos_dict
        )
        args_dict = self.input_dict['dos']['args']
        args_type = self.input_dict['dos']['args_type']
        input_dos_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_dos_dict
        )

        #Write.
        self.input_wfndos: str = QePwInputFile(input_wfndos_dict, self.input_dict).get_input_str()
        self.input_dos: str = QePwInputFile.write_general(input_dos_dict)

    def set_jobs_str(self):
        self.job_wfndos = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x < wfndos.in &> wfndos.in.out 
'''     
        self.job_dos = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}dos.x -pd .true. < dos.in &> dos.in.out 
'''
        self.jobs = [
            './job_wfndos.sh',
            './job_dos.sh'
        ]

    def create(self):
        write_str_2_f('wfndos.in', self.input_wfndos)
        write_str_2_f('job_wfndos.sh', self.job_wfndos)
        write_str_2_f('dos.in', self.input_dos)
        write_str_2_f('job_dos.sh', self.job_dos)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wfndos.in*',
            'dos.in*',
            'struct_dos.dat',
            'job_dos.sh',
            'job_wfndos.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wfndos.in')
        os.system('rm -rf dos.in')
        os.system('rm -rf job_dos.sh')
        os.system('rm -rf job_wfndos.sh')
        
        os.system('rm -rf struct_dos.dat')
        os.system('rm -rf dos.in.out')
        os.system('rm -rf wfndos.in.out')
#endregion
