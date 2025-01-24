#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EpwJob:
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
        if isinstance(self.input_dict['epw']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['epw']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['epw']['job_info'])
    
    def set_inputs_str(self):
        # Base.
        input_epw_dict = {
            'namelists': {
                'inputepw': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'nk1': self.input_dict['dfpt']['qdim'][0],
                    'nk2': self.input_dict['dfpt']['qdim'][1],
                    'nk3': self.input_dict['dfpt']['qdim'][2],
                    'nq1': self.input_dict['dfpt']['qdim'][0],
                    'nq2': self.input_dict['dfpt']['qdim'][1],
                    'nq3': self.input_dict['dfpt']['qdim'][2],
                    'nkf1': self.input_dict['dfpt']['qdim'][0],
                    'nkf2': self.input_dict['dfpt']['qdim'][1],
                    'nkf3': self.input_dict['dfpt']['qdim'][2],
                    'nqf1': self.input_dict['dfpt']['qdim'][0],
                    'nqf2': self.input_dict['dfpt']['qdim'][1],
                    'nqf3': self.input_dict['dfpt']['qdim'][2],
                    'nbndsub': self.input_dict['abs']['num_val_bands'] + self.input_dict['abs']['num_cond_bands'],
                    'dvscf_dir': "'./save'",
                    'elph': '.true.',
                    'epbwrite': '.true.',
                    'epbread': '.false.',
                    'prtgkk': '.false.',
                    'wannierize': '.true.',
                    'auto_projections': '.true.',
                    'scdm_proj': '.true.',
                    # 'temps': '300.0',
                    # 'verbosity': 1,
                }
            }
        }

        # Additions.
        exclude_bands_str = self.input.epw.get_skipped_bands_str()
        if exclude_bands_str is not None: input_epw_dict['namelists']['inputepw']['bands_skipped'] = exclude_bands_str
        args_dict = self.input_dict['epw']['args']
        args_type = self.input_dict['epw']['args_type']
        input_epw_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_epw_dict
        )

        # Write.
        self.input_epw: str = QePwInputFile.write_general(input_epw_dict)

    def set_jobs_str(self):
        self.job_epw = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}{self.input_dict['epw']['exec_loc']} {self.scheduler.get_sched_mpi_infix(self.job_info)} < epw.in  &> epw.in.out 
cp ./wfn.xml ./save/wfn.xml
cp ./input_qpt_update.pkl ./save/input.pkl
cp ./tmp/*epb* ./save/
'''
        
        self.jobs = [
            './job_epw.sh',
        ]

    def create(self):
        write_str_2_f('epw.in', self.input_epw)
        write_str_2_f('job_epw.sh', self.job_epw)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'epw.in*',
            'save',
            'job_epw.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf epw.in')
        os.system('rm -rf job_epw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf ./struct*')
        os.system('rm -rf ./decay*')
        os.system('rm -rf ./struct_elph*')
        os.system('rm -rf EPW.bib')
        os.system('rm -rf epwdata.fmt')
        os.system('rm -rf selecq.fmt')
        os.system('rm -rf vmedata.fmt')
        os.system('rm -rf crystal.fmt')
        os.system('rm -rf epw.in.out')
        os.system('rm -rf epw.out')
#endregion
