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
class Phdos:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_q2r_dos = \
f'''&INPUT
zasr='crystal'
fildyn='struct.dyn'
flfrc='struct.fc'
/
'''
        
        self.job_q2r_dos = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.phdos.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.phdos.job_desc)}q2r.x < q2r_dos.in &> q2r_dos.in.out 
'''
        
        self.input_matdyn_dos = \
f'''&INPUT
asr='crystal'
flfrc='struct.fc'
flfrq='struct.phdos.freq'
flvec='struct.phdos.modes'
dos=.true.
fldos='struct.phdos'
nk1={int(self.input.phdos.qdim[0])}
nk2={int(self.input.phdos.qdim[1])}
nk3={int(self.input.phdos.qdim[2])}
/
'''
        
        self.job_matdyn_dos = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.phbands.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.phbands.job_desc)}matdyn.x < matdyn_dos.in &> matdyn_dos.in.out 
'''

        self.jobs = [
            'job_q2r_dos.sh',
            'job_matdyn_dos.sh',
        ]

    def create(self):
        write_str_2_f('q2r_dos.in', self.input_q2r_dos)
        write_str_2_f('job_q2r_dos.sh', self.job_q2r_dos)
        write_str_2_f('matdyn_dos.in', self.input_matdyn_dos)
        write_str_2_f('job_matdyn_dos.sh', self.job_matdyn_dos)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_q2r_dos.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_matdyn_dos.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'q2r_dos.in*',
            'job_q2r_dos.sh',
            'matdyn_dos.in*',
            'job_matdyn_dos.sh',
            'struct.dyn*',
            'struct.fc',
            'struct.freq*',
            'struct.phdos*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf q2r_dos.in')
        os.system('rm -rf q2r_dos.in.out')
        os.system('rm -rf job_q2r_dos.sh')
        
        os.system('rm -rf matdyn_dos.in')
        os.system('rm -rf matdyn_dos.in.out')
        os.system('rm -rf job_matdyn_dos.sh')
        
        os.system('rm -rf struct.dyn*')
        os.system('rm -rf struct.fc')
        os.system('rm -rf struct.freq')
        os.system('rm -rf struct.freq.gp')
        os.system('rm -rf struct.phdos')
        os.system('rm -rf struct.phdos.freq')
        os.system('rm -rf struct.phdos.freq.gp')
        os.system('rm -rf struct.phdos.modes')
#endregion
