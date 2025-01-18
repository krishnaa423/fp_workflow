#region: Modules.
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

    def set_update_files(self):
        self.input_update_qpts: str = '''

#region modules
import xml.etree.ElementTree as ET
from io import StringIO
from fp.io.pkl import load_obj, save_obj
from fp.inputs.input_main import Input
from fp.jobs.dfpt import DfptJob
from fp.jobs.epw import EpwJob
from fp.schedulers.scheduler import JobProcDesc
from fp.flows.flow_manage import FlowManage
#endregions

#region variables
#endregions

#region functions
def main():
    scf_updater = ScfUpdater()
    scf_updater.update()
#endregions

#region classes
class ScfUpdater:
    def __init__(self):
        self.nk: int = None
        self.input: Input = None
        self.flowmanage: FlowManage = None
        
    def read_nk(self):
        with open('./scf.xml', 'r') as f: 
            root = ET.parse(f).getroot()
        self.nk = int(root.findall('.//nks')[0].text)
    
    def read_input_and_flow(self):
        self.input = load_obj('./input.pkl')
        self.flowmanage = load_obj('./flowmanage.pkl')

    def update(self):
        # Read.
        self.read_nk()
        self.read_input_and_flow()

        update_files: list = self.input.input_dict['scf']['update_files']

        for list_step in self.flowmanage.list_of_steps:
            # dfpt.
            if 'job_dfpt.sh' in update_files and isinstance(list_step, DfptJob):
                # Update.
                dfpt: DfptJob = DfptJob(self.input)
                job_info: JobProcDesc = dfpt.job_info
                if job_info.ni is None:
                    job_info.ni = self.nk
                    self.input.input_dict['dfpt']['job_info'] = {
                        'nodes': job_info.nodes,
                        'ntasks': job_info.ntasks,
                        'time': job_info.time,
                        'ni': job_info.ni,
                        'nk': job_info.nk,
                    }
            
                    # Write.
                    dfpt: DfptJob = DfptJob(self.input)
                    dfpt.create()
                    save_obj(self.input, 'input_qpt_update.pkl')
        
#endregions

#region main
main()
#endregions
'''

    def set_inputs_str(self):
        self.set_update_files()
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

    def get_update_job_cmd(self) -> str:
        output = ''
        if 'job_dfpt.sh' in self.input.input_dict['scf']['update_files']:
            output = 'python3 update_numqpts_from_scf.py\n'

        return output

    def set_jobs_str(self):
        update_output_str: str = self.get_update_job_cmd()
        self.job_scf: str = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < scf.in &> scf.in.out

cp ./tmp/struct.save/data-file-schema.xml ./scf.xml
{update_output_str}
'''
    
        self.jobs = [
            './job_scf.sh',
        ]

    def create_update_files(self):
        if self.input.input_dict.get('scf', {}).get('update_files') is not None:
            write_str_2_f('update_numqpts_from_scf.py', self.input_update_qpts)

    def create(self):
        self.create_update_files()
        write_str_2_f('scf.in', self.input_scf)
        write_str_2_f('job_scf.sh', self.job_scf)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
        os.system('rm -rf update_numqpts_from_scf.py')
        os.system('rm -rf input_qpt_update.pkl')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf scf.in.out')
        os.system('rm -rf scf.xml')
        os.system('rm -rf pseudos')

#endregion
