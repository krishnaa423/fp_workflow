#region modules
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
from fp.inputs.bgw import BgwInputFile
from fp.flows.fullgridflow import FullGridFlow
from fp.flows.flow_manage import FlowManage
from fp.jobs.scf import ScfJob
from fp.jobs.dfpt import DfptJob
from fp.jobs.phbands import PhbandsJob
from fp.jobs.phdos import PhdosJob
from fp.jobs.phmodes import PhmodesJob
from fp.jobs.dftelbands import DftelbandsJob
from fp.jobs.wfngeneral import WfnJob, WfnqJob, WfnfiJob, WfnqfiJob
from fp.jobs.epsilon import EpsilonJob
from fp.jobs.sigma import SigmaJob
from fp.jobs.inteqp import InteqpJob
from fp.jobs.kernel import KernelJob
from fp.jobs.abs import AbsorptionJob
import copy 
#endregions

#region variables
jobclass_map = {
    'ScfJob': ScfJob,
    'DfptJob': DfptJob,
    'PhbandsJob': PhbandsJob,
    'PhmodesJob': PhmodesJob,
    'PhdosJob': PhdosJob,
    'DftelbandsJob': DftelbandsJob,
    'WfnJob': WfnJob,
    'WfnqJob': WfnqJob,
    'WfnfiJob': WfnfiJob,
    'WfnqfiJob': WfnqfiJob,
    'EpsilonJob': EpsilonJob,
    'SigmaJob': SigmaJob,
    'KernelJob': KernelJob,
    'AbsorptionJob': AbsorptionJob,
}
#endregions

#region functions
#endregions

#region classes
class ConvergenceJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)

        self.conv_dict: dict = self.input_dict.get('convergence')


        # For now, will not populate this. Convergence tasks need to be done seperately.
        self.jobs: list = []

        # Shared vars.
        self.folders: dict = {
            'dft': [],
            'dfpt': [],
            'gw': [],
            'bse': [],
        }

        self.populate_job_list()

    def populate_job_list(self):
        conv_types = [
            'dft',
            'dfpt',
            'gw',
            'bse',
        ]

        for conv_type in conv_types:
            if self.conv_dict.get(conv_type) is not None:
                for conv_dict_list_item  in self.conv_dict.get(conv_type):
                    # Get directory.
                    dir_suffix = ''
                    for key, value in conv_dict_list_item.items():
                        if key!= 'reuse':
                            for value_key, value_value in value.items():
                                if isinstance(value_value, list):
                                    value_value = 'x'.join([str(item) for item in value_value])
                                dir_suffix += f'{key}@{value_key}~{value_value}#'
                    dest_directory = f'./convergence/{conv_type}/{conv_type}={dir_suffix}'
                    
                    # Add directory. 
                    if self.conv_dict.get('add_to_job_list') is not None:
                        if self.conv_dict['add_to_job_list']:
                            self.jobs.append([
                                dest_directory,
                                './runall.sh',
                            ])

    def do_reuse(self, current_folder: str, reuse_dict: dict):
        # Can have linking to many folders.
        for copy_files_dict in reuse_dict['copy_files']:
            folder_type = copy_files_dict['folder']
            subfolder = self.folders[folder_type][copy_files_dict['idx']]
            source_folder = f'../../{folder_type}/{subfolder}'

            # Copy files from this folder.
            for file in copy_files_dict['files']:
                source_file = os.path.join(source_folder, file)
                dest_file = os.path.join(current_folder, file)
                os.system(f'ln -sf {source_file} {dest_file}')

    def remove_jobclasses_from_list(self, jobnamelist_to_remove: list, jobclasslist_to_remove_from: list) -> list:
        output_list = copy.deepcopy(jobclasslist_to_remove_from)
        
        for jobname_item in jobnamelist_to_remove:
            output_list.remove(jobclass_map[jobname_item])

        return output_list

    def create_generic(
        self,
        conv_type: str,
        list_of_steps: list,
        start_idx: str,
        stop_idx: str,
    ):
        os.system(f'mkdir -p ./convergence/{conv_type}')

        # Iterate.
        if self.conv_dict.get(conv_type) is not None:
            for conv_dict_list_item  in self.conv_dict.get(conv_type):
                # Prefix.
                input_dict_copy = copy.deepcopy(self.input_dict)      # So that we don't change it for other functions.
                dir_suffix = ''
                for key, value in conv_dict_list_item.items():
                    if key!= 'reuse':
                        for value_key, value_value in value.items():
                            if isinstance(value_value, list):
                                value_value = 'x'.join([str(item) for item in value_value])
                            dir_suffix += f'{key}@{value_key}~{value_value}#'
                dest_directory = f'./convergence/{conv_type}/{conv_type}={dir_suffix}'
                self.folders[conv_type].append(f'{conv_type}={dir_suffix}')
                os.system(f'mkdir -p {dest_directory}')

                # Generate.
                #Link files.
                filtered_steps = None
                if conv_dict_list_item.get('reuse') is not None:
                    self.do_reuse(dest_directory, conv_dict_list_item['reuse'])
                    filtered_steps = self.remove_jobclasses_from_list(
                        conv_dict_list_item['reuse']['turn_off_jobs'],
                        list_of_steps
                    )
                    conv_dict_list_item.pop('reuse', None)
                #Job files.  
                for key, value in conv_dict_list_item.items():
                    input_dict_copy[key].update(value)
                FullGridFlow.create_from_list(
                    source_input_dict=input_dict_copy,
                    list_of_steps=filtered_steps if filtered_steps is not None else list_of_steps,
                    start_idx=start_idx if filtered_steps is None else None,
                    stop_idx=stop_idx if filtered_steps is None else None,
                    dest_directory=dest_directory,
                    copy_additional=self.input_dict.get('convergence', {}).get('copy_additional')
                )


    def create_dft(self):
        self.create_generic(
            conv_type='dft',
            list_of_steps=[
                ScfJob,
                DftelbandsJob,
            ],
            start_idx=0,
            stop_idx=1,
        )

    def create_dfpt(self):
        self.create_generic(
            conv_type='dfpt',
            list_of_steps=[
                ScfJob,
                DfptJob,
                PhbandsJob,
                PhdosJob,
                PhmodesJob,
            ],
            start_idx=0,
            stop_idx=3,
        )

    def create_gw(self):
        self.create_generic(
            conv_type='gw',
            list_of_steps=[
                ScfJob,
                DftelbandsJob,
                WfnJob,
                WfnqJob,
                EpsilonJob,
                SigmaJob,
                InteqpJob,
            ],
            start_idx=0,
            stop_idx=10,
        )

    def create_bse(self):
        self.create_generic(
            conv_type='bse',
            list_of_steps=[
                ScfJob,
                WfnJob,
                WfnqJob,
                WfnfiJob,
                WfnqfiJob,
                EpsilonJob,
                SigmaJob,
                KernelJob,
                AbsorptionJob,
            ],
            start_idx=0,
            stop_idx=13,
        )

    def create(self):
        os.system('mkdir -p ./convergence')
        self.create_dft()
        self.create_dfpt()
        self.create_gw()
        self.create_bse()
        
    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'epsilon.inp*',
            'job_epsilon.sh',
            'epsmat.h5',
            'eps0mat.h5',
            'epsilon.log',
            'chi_converge.dat',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf convergence')

#endregions
