#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
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
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['dfpt']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['dfpt']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['dfpt']['job_info'])

        # Set to the max multiple for ntasks w.r.t. ni. 
        # TODO: nk. Currently only image parallelization through -ni flag. 
        if self.job_info.ni is not None:
            ntasks = self.job_info.ntasks
            ni = self.job_info.ni
            if ntasks < ni: 
                ntasks = ni
                self.job_info.ntasks = ntasks 
            elif ntasks % ni != 0:
                updated_ntasks = ( ntasks // ni ) * ni 
                # Set job_info tasks.
                self.job_info.ntasks = updated_ntasks

        self.job_recover_info = JobProcDesc(
            nodes=self.job_info.nodes,
            ntasks=self.job_info.ntasks,
            time=self.job_info.time,
            ni=self.job_info.ni,
            nk=self.job_info.nk,
        )
        if self.job_recover_info.ni is not None:
            self.job_recover_info.ntasks /= self.job_recover_info.ni
            self.job_recover_info.ntasks = int(self.job_recover_info.ntasks)
            self.job_recover_info.ni = None

    def set_inputs_str(self):
        # Base. 
        input_dfpt_dict = {
            'namelists': {
                'inputph': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'ldisp': '.true.',
                    'nq1': self.input_dict['dfpt']['qdim'][0],
                    'nq2': self.input_dict['dfpt']['qdim'][1],
                    'nq3': self.input_dict['dfpt']['qdim'][2],
                    'fildyn': "'struct.dyn'",
                    'tr2_ph': self.input_dict['dfpt']['conv_threshold'],
                    'fildvscf': "'dvscf'",
                }
            }
        }

        # Additions.
        #override or extra. 
        args_dict = self.input_dict['dfpt']['args']
        args_type = self.input_dict['dfpt']['args_type']
        input_dfpt_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_dfpt_dict
        )

        input_dfpt_recover_dict = input_dfpt_dict.copy()
        input_dfpt_recover_dict['namelists']['inputph']['recover'] = '.true.'

        # Get inputs.
        self.input_dfpt: str = QePwInputFile.write_general(input_dfpt_dict)
        self.input_dfpt_recover: str = QePwInputFile.write_general(input_dfpt_recover_dict)

    def set_jobs_str(self):
        self.jobs = [
            './job_dfpt.sh',
        ]

        self.job_dfpt = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}ph.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < dfpt.in &> dfpt.in.out  
{self.scheduler.get_sched_mpi_prefix(self.job_recover_info)}ph.x {self.scheduler.get_sched_mpi_infix(self.job_recover_info)} < dfpt_recover.in &> dfpt_recover.in.out  

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
        os.system('rm -rf dfpt*.in')
        os.system('rm -rf dfpt*.in.out')
        os.system('rm -rf create_save.py')
        os.system('rm -rf job_dfpt.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf out*')
        os.system('rm -rf ./save')
        os.system('rm -rf struct.dyn*')

#endregion
