#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
#endregion

#region: Variables.
update_from_relax_file_content = '''#!/usr/bin/env python3

#region modules
from fp.inputs.qepw import QePwInputFile
from fp.io.strings import write_str_2_f
from fp.io.pkl import load_obj
from fp.inputs.atoms import AtomsInput
from typing import List
import os
from fp.flows.flow_manage import FlowManage
from fp.jobs.phonopy import PhonopyJob
from fp.jobs.wannier import WannierJob
from fp.inputs.input_main import Input
#endregions

#region variables
#endregions

#region functions
def main():
    relax_updater = RelaxUpdater()
    relax_updater.update()
#endregions

#region classes
class RelaxUpdater:
    def __init__(self):
        # Read from input_dict.
        self.input_dict: dict = load_obj('./input_dict.pkl')
        self.flowmanage: FlowManage = load_obj('./flowmanage.pkl')

        # Vars across methods.
        self.input: Input = None
        self.update_input()
        self.atoms_input = self.input.atoms

        # Files to update.
        self.files: List[str] = self.input_dict.get('relax', {}).get('update_files')

    def update_input(self):
        self.input_dict['atoms']['read']['cell']['file'] = 'relaxed_cell_parameters.txt'
        self.input_dict['atoms']['read']['positions']['file'] = 'relaxed_atomic_positions.txt'
        self.input: Input = Input.from_dict(self.input_dict)
    
    def get_pw(self, filename: str) -> dict:
        with open(filename, 'r') as f: file_content = f.read()
        
        return QePwInputFile.read_pw(file_content)

    def update(self):
        if self.files is None: return

        for file in self.files:
            if os.path.isfile(file):
                pw_dict: dict = self.get_pw(file)
                pw_dict['blocks']['atomic_positions'] = self.atoms_input.get_qe_scf_atomic_positions()
                pw_dict['blocks']['cell_parameters'] = self.atoms_input.get_qe_scf_cell()
                pw_writer = QePwInputFile(pw_dict, self.input_dict)
                write_str_2_f(file, pw_writer.get_input_str())

        
        # Update phonopy.
        for list_step in self.flowmanage.list_of_steps:
            if isinstance(list_step, PhonopyJob):
                PhonopyJob(self.input).create()
            if isinstance(list_step, WannierJob):
                WannierJob(self.input).create()
    
#endregions

#region main
main()
#endregions

'''
#endregion

#region: Functions.
#endregion

#region: Classes.
class RelaxJob:
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
        if isinstance(self.input_dict['relax']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['relax']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['relax']['job_info'])

    def set_inputs_str(self):
        #Base. 
        input_relax_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': f'{self.input.relax.calc_str()}',
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
        #occupations.
        if 'cdft' in self.input_dict['relax']['type']: 
            input_relax_dict['namelists']['system']['occupations'] = 'from_input'
            input_relax_dict['blocks']['occupations'] = self.input.relax.get_occupations()
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_relax_dict['namelists']['system']['noncolin'] = True
            input_relax_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['relax']['args']
        args_type = self.input_dict['relax']['args_type']
        input_relax_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_relax_dict
        )

        # Get string. 
        self.input_relax: str = QePwInputFile(input_relax_dict, self.input_dict).get_input_str()

    def set_jobs_str(self):
        cell_unit_str = self.input_dict['atoms']['write_cell_units']
        pos_unit_str = self.input_dict['atoms']['write_position_units']
        
        # Cell parameters. 
        save_final_cell_parameters_str = \
"awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /CELL_PARAMETERS/ {cell_flag=1; next} /End final coordinates/ {end_flag=0} end_flag && cell_flag {print; if (length==0) cell_flag=0 }' relax.in.out > relaxed_cell_parameters.txt"
        save_final_cell_parameters_str += f'\necho "{cell_unit_str}" | cat - relaxed_cell_parameters.txt > temp && mv temp relaxed_cell_parameters.txt'

        # Atomic positions.
        save_final_atomic_positions_str = \
"awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /ATOMIC_POSITIONS/ {pos_flag=1; next} /End final coordinates/ {end_flag=0}  end_flag && pos_flag { print $1, $2, $3, $4 }' relax.in.out > relaxed_atomic_positions.txt"
        save_final_atomic_positions_str += f'\necho "{pos_unit_str}" | cat - relaxed_atomic_positions.txt > temp && mv temp relaxed_atomic_positions.txt'

        self.job_relax: str = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < relax.in &> relax.in.out

cp ./tmp/struct.save/data-file-schema.xml ./relax.xml

# Copy the end atomic positions and cell parameters (if vc-relax).
{save_final_cell_parameters_str}
{save_final_atomic_positions_str}

# Update from relax.
./update_from_relax.sh
'''
    
        self.jobs = [
            './job_relax.sh',
        ]

    def create(self):
        write_str_2_f('relax.in', self.input_relax)
        write_str_2_f('job_relax.sh', self.job_relax)
        write_str_2_f('update_from_relax.sh', update_from_relax_file_content)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'relax.in',
            'job_relax.sh',
            'relax.in.out',
            'relax.xml',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf relax.in')
        os.system('rm -rf job_relax.sh')
        os.system('rm -rf update_from_relax.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf relax.in.out')
        os.system('rm -rf relax.xml')

#endregion
