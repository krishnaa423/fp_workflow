#region: Modules.
from fp.inputs.input_main import *
from fp.io.strings import *
from fp.flows.run import *
from pkg_resources import resource_filename
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DfptJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
    
        self.jobs = [
            'job_dfpt.sh',
        ]

        self.input_dfpt = \
f'''&INPUTPH
outdir='./tmp'
prefix='struct'
ldisp=.true.
nq1={self.input.dfpt.qgrid[0]}
nq2={self.input.dfpt.qgrid[1]}
nq3={self.input.dfpt.qgrid[2]}
tr2_ph={self.input.dfpt.conv_threshold}
fildyn='struct.dyn'
fildvscf='dvscf'
{self.input.dfpt.extra_args if self.input.dfpt.extra_args is not None else ""}
/
'''
        self.input_dfpt_recover = \
f'''&INPUTPH
outdir='./tmp'
prefix='struct'
ldisp=.true.
nq1={self.input.dfpt.qgrid[0]}
nq2={self.input.dfpt.qgrid[1]}
nq3={self.input.dfpt.qgrid[2]}
tr2_ph={self.input.dfpt.conv_threshold}
fildyn='struct.dyn'
fildvscf='dvscf'
recover=.true.
{self.input.dfpt.extra_args if self.input.dfpt.extra_args is not None else ""}
/
'''
        dfpt_recover_job_desc = JobProcDesc(
            nodes=self.input.dfpt.job_desc.nodes,
            ntasks=self.input.dfpt.job_desc.ntasks,
            time=self.input.dfpt.job_desc.time,
            ni=self.input.dfpt.job_desc.ni,
            nk=self.input.dfpt.job_desc.nk,
        )
        if dfpt_recover_job_desc.ni:
            dfpt_recover_job_desc.ntasks /= dfpt_recover_job_desc.ni
            dfpt_recover_job_desc.ntasks = int(dfpt_recover_job_desc.ntasks)

        self.job_dfpt = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dfpt.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dfpt.job_desc)}ph.x {self.input.scheduler.get_sched_mpi_infix(self.input.dfpt.job_desc)} < dfpt.in &> dfpt.in.out  
{self.input.scheduler.get_sched_mpi_prefix(dfpt_recover_job_desc)}ph.x {self.input.scheduler.get_sched_mpi_infix(dfpt_recover_job_desc, add_ni_if_present=False)} < dfpt_recover.in &> dfpt_recover.in.out  

python3 ./create_save.py
'''

    def copy_createsave_file(self):
        pkg_dir = resource_filename('fp', '')
        src_path = pkg_dir + '/jobs/create_save.py'
        dst_path = './create_save.py'

        os.system(f'cp {src_path} {dst_path}')

    def create(self):
        self.copy_createsave_file()

        write_str_2_f('dfpt.in', self.input_dfpt)
        write_str_2_f('dfpt_recover.in', self.input_dfpt_recover)
        write_str_2_f('job_dfpt.sh', self.job_dfpt)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'dfpt*.in',
            'dfpt*.in.out',
            'job_dfpt.sh',
            'save',
            'struct.dyn*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf dfpt_start.in')
        os.system('rm -rf dfpt_start.in.out')
        os.system('rm -rf dfpt_end.in')
        os.system('rm -rf dfpt_end.in.out')
        os.system('rm -rf dfpt*.in')
        os.system('rm -rf dfpt*.in.out')
        os.system('rm -rf create_save.py')
        os.system('rm -rf job_dfpt.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf out*')
        os.system('rm -rf ./save')
        os.system('rm -rf struct.dyn*')
#endregion
