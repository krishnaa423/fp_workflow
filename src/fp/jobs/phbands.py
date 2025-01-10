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
class PhbandsJob:
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
        if isinstance(self.input_dict['phbands']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['phbands']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['phbands']['job_info'])

    def set_inputs_str(self):
        # Base. 
        input_q2r_bands_dict = {
            'namelists': {
                'input': {
                    'zasr': "'crystal'",
                    'fildyn': "'struct.dyn'",
                    'flfrc': "'struct.fc'",
                }
            }
        }
        input_matdyn_bands_dict = {
            'namelists': {
                'input': {
                    'asr': "'crystal'",
                    'flfrc': "'struct.fc'",
                    'flfrq': "'struct.freq'",
                    'flvec': "'struct.modes'",
                    'q_in_band_form': '.true.',
                    'q_in_cryst_coord': '.true.',
                }
            }
        }

        # Additions.
        #q2r.
        args_dict = self.input_dict['phbands']['q2r_args']
        args_type = self.input_dict['phbands']['q2r_args_type']
        input_q2r_bands_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_q2r_bands_dict
        )
        #matdyn.
        args_dict = self.input_dict['phbands']['matdyn_args']
        args_type = self.input_dict['phbands']['matdyn_args_type']
        input_matdyn_bands_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_matdyn_bands_dict
        )

        # Get inputs.
        self.input_q2r_bands: str = QePwInputFile.write_general(input_q2r_bands_dict)
        self.input_matdyn_bands: str = QePwInputFile.write_general(input_matdyn_bands_dict)
        self.input_matdyn_bands += f'{self.input.phbands.get_kpath_str()}'

    def set_jobs_str(self):
        
        self.job_q2r_bands = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}q2r.x < q2r_bands.in &> q2r_bands.in.out 
'''
        
        self.job_matdyn_bands = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}matdyn.x < matdyn_bands.in &> matdyn_bands.in.out 
'''

        self.jobs = [
            './job_q2r_bands.sh',
            './job_matdyn_bands.sh',
        ]

    def create(self):
        write_str_2_f('q2r_bands.in', self.input_q2r_bands)
        write_str_2_f('job_q2r_bands.sh', self.job_q2r_bands)
        write_str_2_f('matdyn_bands.in', self.input_matdyn_bands)
        write_str_2_f('job_matdyn_bands.sh', self.job_matdyn_bands)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'q2r_bands.in*',
            'job_q2r_bands.sh',
            'matdyn_bands.in*',
            'job_matdyn_bands.sh',
            'struct.dyn*',
            'struct.fc',
            'struct.freq',
            'struct.modes',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf q2r_bands.in')
        os.system('rm -rf q2r_bands.in.out')
        os.system('rm -rf job_q2r_bands.sh')
        
        os.system('rm -rf matdyn_bands.in')
        os.system('rm -rf matdyn_bands.in.out')
        os.system('rm -rf job_matdyn_bands.sh')
        
        os.system('rm -rf struct.dyn*')
        os.system('rm -rf struct.fc')
        os.system('rm -rf struct.freq')
        os.system('rm -rf struct.modes')
#endregion
