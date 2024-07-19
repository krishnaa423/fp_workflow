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
class Dfpt:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
    
        self.jobs = [
            'job_dfpt.sh',
        ]

    def copy_createsave_file(self):
        pkg_dir = resource_filename('fp', '')
        src_path = pkg_dir + '/calcs/create_save.py'
        dst_path = './create_save.py'

        os.system(f'cp {src_path} {dst_path}')

    def create_one_job(self):
        input_dfpt = \
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
/
'''

        job_dfpt = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dfpt.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dfpt.job_desc)}ph.x {self.input.scheduler.get_sched_mpi_infix(self.input.dfpt.job_desc)} < dfpt.in &> dfpt.in.out  

python3 ./create_save.py
'''

        write_str_2_f('dfpt.in', input_dfpt)
        write_str_2_f('job_dfpt.sh', job_dfpt)

    def get_input_njob_file(self, start_irr, last_irr):
        output = \
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
start_irr={start_irr}
last_irr={last_irr}
'''
        return output 

    def get_inputs_middle(self):
        output_list = []
        for job_idx in range(self.input.dfpt.njobs):
            output_list.append(self.get_input_njob_file(*self.input.dfpt.get_irr_lim(job_idx)))

        return output_list

    def get_job_dfpt(self):
        njobs = self.input.dfpt.njobs

        script_middle_block = ''
        for job_idx in range(njobs):
            total_tasks = self.input.dfpt.job_desc.ntasks 
            job_tasks = total_tasks // njobs if job_idx != njobs-1 else total_tasks - (total_tasks // njobs)*(njobs-1)
            split_job_desc = JobProcDesc(
                time=self.input.dfpt.job_desc.time,
                nodes=self.input.dfpt.job_desc.nodes,
                ntasks=job_tasks,
                nk=self.input.dfpt.job_desc.nk,
                ni=self.input.dfpt.job_desc.ni,
            )

            script_middle_block += f'{self.input.scheduler.get_sched_mpi_prefix(split_job_desc)}ph.x {self.input.scheduler.get_sched_mpi_infix(split_job_desc)} < dfpt_{job_idx}.in &> dfpt_{job_idx}.in.out &\n'

        output = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dfpt.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dfpt.job_desc)}ph.x {self.input.scheduler.get_sched_mpi_infix(self.input.dfpt.job_desc)} < dfpt_start.in &> dfpt_start.in.out  
{script_middle_block}
wait
{self.input.scheduler.get_sched_mpi_prefix(self.input.dfpt.job_desc)}ph.x {self.input.scheduler.get_sched_mpi_infix(self.input.dfpt.job_desc)} < dfpt_end.in &> dfpt_end.in.out  
python3 create_save.py
'''
        
        return output

    def create_njobs(self):
        input_start = \
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
start_irr=0
last_irr=0
/
'''
        
        inputs_middle = self.get_inputs_middle()

        input_end = \
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
/
'''
        
        job_dfpt = self.get_job_dfpt()


        write_str_2_f('dfpt_start.in', input_start)
        for job_idx, input_middle in zip(range(self.input.dfpt.njobs), inputs_middle):
            write_str_2_f(f'dfpt_{int(job_idx)}.in', input_middle)
        write_str_2_f('dfpt_end.in', input_end)
        write_str_2_f('job_dfpt.sh', job_dfpt)

    def create(self):
        if self.input.dfpt.njobs:
            self.create_njobs()
        else:
            self.create_one_job()

        self.copy_createsave_file()

    def run(self, total_time):
        total_time = run_and_wait_command('./job_dfpt.sh', self.input, total_time)

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
        os.system('rm -rf ./save')
        os.system('rm -rf struct.dyn*')
        

        # if self.input.dfpt.njobs: 
        #     for job_idx in range(self.input.dfpt.njobs):
        #         os.system(f'rm -rf dfpt_{job_idx}.in')
        #         os.system(f'rm -rf dfpt_{job_idx}.in.out')
#endregion
