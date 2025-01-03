#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.jobs.qepw import QePwInputFile, IbravType
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class ScfJob:
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
        if isinstance(self.input_dict['scf']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['scf']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['scf']['job_info'])

    def set_inputs_str(self):
        #Base. 
        input_scf_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'scf',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
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
                'kpoints': self.input_dict['scf']['kdim'],
            },
            'kpoints_type': 'automatic',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_scf_dict['namelists']['system']['noncolin'] = True
            input_scf_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['scf']['args']
        args_type = self.input_dict['scf']['args_type']
        input_scf_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_scf_dict
        )

        # Get string. 
        self.input_scf: str = QePwInputFile(input_scf_dict, self.input_dict).get_input_str()

    def set_jobs_str(self):
        self.job_scf: str = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < scf.in &> scf.in.out

cp ./tmp/struct.save/data-file-schema.xml ./scf.xml
'''
    
        self.jobs = [
            'job_scf.sh',
        ]

    def create(self):
        write_str_2_f('scf.in', self.input_scf)
        write_str_2_f('job_scf.sh', self.job_scf)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_scf.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'scf.in',
            'job_scf.sh',
            'tmp',
            'scf.xml',
            '*.pkl',
            '*.xsf',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf scf.in')
        os.system('rm -rf job_scf.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf scf.in.out')
        os.system('rm -rf scf.xml')
        os.system('rm -rf pseudos')
#endregion
