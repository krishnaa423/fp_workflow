#region: Modules.
from fp.inputs import *
from fp.io import *
from fp.flows import *
from pkg_resources import resource_filename
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Phbands:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_q2r_bands = \
f'''&INPUT
zasr='crystal'
fildyn='struct.dyn'
flfrc='struct.fc'
/
'''
        
        self.job_q2r_bands = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.phbands.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.phbands.job_desc)}q2r.x < q2r_bands.in &> q2r_bands.in.out 
'''
        
        self.input_matdyn_bands = \
f'''&INPUT
asr='crystal'
flfrc='struct.fc'
flfrq='struct.freq'
flvec='struct.modes'
q_in_cryst_coord=.true.
/
{self.input.phbands.get_kpath_str()}
'''
        
        self.job_matdyn_bands = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.phbands.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.phbands.job_desc)}matdyn.x < matdyn_bands.in &> matdyn_bands.in.out 
'''

        self.jobs = [
            'job_q2r_bands.sh',
            'job_matdyn_bands.sh',
        ]

    def create(self):
        write_str_2_f('q2r_bands.in', self.input_q2r_bands)
        write_str_2_f('job_q2r_bands.sh', self.job_q2r_bands)
        write_str_2_f('matdyn_bands.in', self.input_matdyn_bands)
        write_str_2_f('job_matdyn_bands.sh', self.job_matdyn_bands)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_q2r_bands.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_matdyn_bands.sh', self.input, total_time)

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
