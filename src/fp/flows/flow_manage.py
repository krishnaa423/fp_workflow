#region: Modules.
from ase import Atoms 
from importlib.util import find_spec
import os 
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class FlowManage:
    def __init__(
        self,
        list_of_steps,
    ):
        self.list_of_steps: list = list_of_steps

    @staticmethod
    def create_pseudos(atoms: Atoms, is_fr: bool=False, xc_type: str='pbe', pseudos_dict: dict=None):
        # Flags. 
        sr_fr_str = 'fr' if is_fr else 'sr'

        # Directories. 
        pkg_dir = os.path.dirname(find_spec('fp').origin)
        pseudo_dir = pkg_dir + f'/data/pseudos/qe/{sr_fr_str}_{xc_type}'
        os.system('mkdir -p ./pseudos')

        # Copy pseudos.
        symbols = atoms.get_chemical_symbols()
        for sym in symbols:
            source_file = pseudo_dir + f'/{sym}.upf'
            dest_file = './pseudos' + f'/{sym}.upf'
            os.system(f'cp {source_file} {dest_file}')

        # Override pseudos.
        if pseudos_dict is not None:
            if pseudos_dict.get('override') is not None:
                for item in pseudos_dict['override']:
                    os.system(f'cp {item["source_fileloc"]} ./pseudos/{item["dest_filename"]}')

        # Generate and add doped pseudopotentials if needed.

    def create_jobs(self):
        assert len(self.list_of_steps)>=1, 'Number of steps/jobs should be greater than 0.'

        for step in self.list_of_steps:
            step.create()

    def run(self, total_time=0.0, start_idx=None, stop_idx=None):
        assert len(self.list_of_steps)>=1 , 'List of steps must have atleast one element.'
        
        total_time: float = total_time

        job_info_string = ''
        job_list= []
        dir_list=  []
        job_idx_counter = 0
        for step in self.list_of_steps:
            for job_info in step.jobs:
                if isinstance(job_info, list): # Has directory and then script.
                    dir_list.append(job_info[0])
                    job_list.append(job_info[1])
                    job_info_string += f'#idx: {job_idx_counter}, dir: {job_info[0]}, script: {job_info[1]}\n'
                else:
                    dir_list.append('./')
                    job_list.append(job_info)
                    job_info_string += f'#idx: {job_idx_counter}, dir: {'./'}, script: {job_info}\n'
                job_idx_counter += 1

        print(f'JOB_INFO:\n{job_info_string}\n\n', flush=True)

        start_idx = start_idx if start_idx is not None else 0
        stop_idx = stop_idx if stop_idx is not None else (job_idx_counter - 1)
        print(f'start_idx: {start_idx}, start_job: {job_list[start_idx]}\n\n\n', flush=True)
        print(f'stop_idx: {stop_idx}, stop_job: {job_list[stop_idx]}\n\n\n', flush=True)

        total_time = 0.0
        for dest_dir, job in zip(dir_list[start_idx:stop_idx+1], job_list[start_idx:stop_idx+1]):
            total_time = run_and_wait_command(f'{job}', self.list_of_steps[0].input, total_time, dest_dir=dest_dir)

        # Write the total workflow run time. 
        print(f'Done whole worflow in {total_time:15.10f} seconds.\n\n', flush=True)

    def save_job_results(self, folder):       
        for step in self.list_of_steps:
            step.save(folder)

    def get_job_all_script(
            self, 
            start_idx=None, 
            stop_idx=None, 
            flowfile='flowmanage.pkl',
        ):
        assert len(self.list_of_steps)>=1, 'There should be atleast one job step.'

        job_info_string = ''
        job_list= []
        dir_list=  []
        job_idx_counter = 0
        for step in self.list_of_steps:
            for job_info in step.jobs:
                if isinstance(job_info, list): # Has directory and then script.
                    dir_list.append(job_info[0])
                    job_list.append(job_info[1])
                    job_info_string += f'#idx: {job_idx_counter}, dir: {job_info[0]}, script: {job_info[1]}\n'
                else:
                    dir_list.append('./')
                    job_list.append(job_info)
                    job_info_string += f'#idx: {job_idx_counter}, dir: {'./'}, script: {job_info}\n'
                job_idx_counter += 1

        output = \
f'''#!/usr/bin/env python3

from fp.flows.flow_manage import FlowManage
import os
from fp.io.pkl import load_obj

# List of jobs:
{job_info_string}

start_idx={start_idx if start_idx is not None else 0}
stop_idx={stop_idx if stop_idx is not None else job_idx_counter-1}

flow: FlowManage = load_obj('{flowfile}')
flow.run(total_time=0, start_idx=start_idx, stop_idx=stop_idx)
'''

        return output 

    def create_job_all_script(self, start_idx=None, stop_idx=None, flowfile_to_read='flowmanage.pkl'):
        write_str_2_f(
            'job_all.sh', 
            self.get_job_all_script(
                start_idx, 
                stop_idx, 
                flowfile_to_read,
            )
        )

        # Write the run script too that wraps around the above script. 
        write_str_2_f(
            'runall.sh',
f'''#!/bin/bash

./job_all.sh &> job_all.out &
'''
        )

    def remove(self, pkl=False, job_all=False, interactive=False, fmt_files=True, xsf=False):
        for step in self.list_of_steps:
            step.remove()

        if pkl: os.system('rm -rf *.pkl')
        if job_all: os.system('rm -rf ./job_all* ./runall.sh')
        if interactive: os.system('rm -rf ./job_interactive.sh')
        if fmt_files: os.system('rm -rf *.fmt')
        if xsf: os.system('rm -rf *.xsf')
                
#endregion