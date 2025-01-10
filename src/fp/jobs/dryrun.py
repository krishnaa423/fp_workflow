#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
import re 
import math 
import subprocess
from fp.schedulers.scheduler import Scheduler, JobProcDesc
import fp.schedulers as schedulers
from fp.inputs.qepw import QePwInputFile, IbravType
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DryrunJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict, runtype='dryrun')
        self.job_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        self.job_info = JobProcDesc(
            nodes=1,
            ntasks=1,
            time='00:20:00',
        )

    def set_inputs_str(self):
        # Base.
        input_dryrun_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': f'md',
                    'nstep': 0,
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
            input_dryrun_dict['namelists']['system']['noncolin'] = True
            input_dryrun_dict['namelists']['system']['lspinorb'] = True
        
        # Get string. 
        self.input_dryrun: str = QePwInputFile(input_dryrun_dict, self.input_dict).get_input_str()

    def set_jobs_str(self):

        # The dryrun jobs will not be in parallel. Just to get some info.
        self.job_dryrun: str = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < dryrun.in &> dryrun.in.out

cp ./tmp/struct.save/data-file-schema.xml ./dryrun.xml
'''

    def create(self):
        write_str_2_f('dryrun.in', self.input_dryrun)
        write_str_2_f('job_dryrun.sh', self.job_dryrun)
        os.system('chmod u+x ./*.sh')

    def run(self, total_time):
        subprocess.run('./job_dryrun.sh')

        return 0.0

    def save(self, folder):
        pass 

    def remove(self):
        os.system('rm -rf dryrun.in')
        os.system('rm -rf job_dryrun.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf dryrun.in.out')
        os.system('rm -rf dryrun.xml')

    def get_max_val(self):
        with open('dryrun.in.out', 'r') as r: txt = r.read()
        pattern  = r'number of electrons\s*=(.*)\n'
        num_of_electrons = int(math.ceil(float(re.findall(pattern, txt)[0])))

        num_bands = int(num_of_electrons/2) if not self.input_dict['scf']['is_spinorbit'] else num_of_electrons
        
        return num_bands
#endregion
