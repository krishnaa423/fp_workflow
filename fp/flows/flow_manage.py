#region: Modules.
from ase import Atoms 
from pkg_resources import resource_filename
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
        is_fr=False,    # Is fully relativistic for the pseudopotentials?
    ):
        self.list_of_steps: list = list_of_steps
        self.is_fr: bool = is_fr

    @staticmethod
    def create_pseudos(atoms: Atoms, is_fr: bool=False):
        sr_fr_str = 'fr' if is_fr else 'sr'

        pkg_dir = resource_filename('fp', '')
        pseudo_dir = pkg_dir + '/data/pseudos/ONCVPSP/sg15'

        os.system('mkdir -p ./ONCVPSP')
        os.system('mkdir -p ./ONCVPSP/sg15')

        symbols = atoms.get_chemical_symbols()

        for sym in symbols:
            source_file = pseudo_dir + f'/{sym}_ONCV_PBE_{sr_fr_str}.upf'
            dest_file = './ONCVPSP/sg15' + f'/{sym}_ONCV_PBE_{sr_fr_str}.upf'
            os.system(f'cp {source_file} {dest_file}')

    def create(self):
        assert len(self.list_of_steps)>=1, 'Number of steps/jobs should be greater than 0.'

        # Add the pseudopotential files.
        self.create_pseudos(self.list_of_steps[0].input.atoms.atoms, self.is_fr)

        for step in self.list_of_steps:
            step.create()

    def run(self, total_time=0.0, start_job=None, stop_job=None):
        assert len(self.list_of_steps)>=1 , 'List of steps must have atleast one element.'
        
        total_time: float = total_time

        job_list = []
        for step in self.list_of_steps:
            for job in step.jobs:
                job_list.append(job)

        print(f'job_list: {job_list}\n\n\n', flush=True)

        start_idx = job_list.index(start_job) if start_job else 0
        stop_idx = job_list.index(stop_job)+1 if stop_job else len(job_list)

        print(f'start_idx: {start_idx}, start_job: {job_list[start_idx]}\n\n\n', flush=True)
        print(f'stop_idx: {stop_idx-1}, stop_job: {job_list[stop_idx-1]}\n\n\n', flush=True)

        total_time = 0.0
        for job in job_list[start_idx:stop_idx]:
            total_time = run_and_wait_command(f'./{job}', self.list_of_steps[0].input, total_time)

        # Write the total workflow run time. 
        print(f'Done whole worflow in {total_time:15.10f} seconds.\n\n', flush=True)

    def save_job_results(self, folder):       
        for step in self.list_of_steps:
            step.save(folder)

    def get_job_all_script(
            self, 
            start_job, 
            stop_job, 
            flowfile='flowmanage.pkl',
        ):
        assert len(self.list_of_steps)>=1, 'There should be atleast one job step.'

        output = \
f'''#!/usr/bin/env python3

from fp.flows import FlowManage
import os
from fp.io import *

start_job='{start_job}'
stop_job='{stop_job}'

flow: FlowManage = load_obj('{flowfile}')
flow.run(total_time=0, start_job=start_job, stop_job=stop_job)
'''

        return output 

    def create_job_all_script(self, filename, start_job, stop_job, flowfile_to_read='flowmanage.pkl'):
        write_str_2_f(
            filename, 
            self.get_job_all_script(
                start_job, 
                stop_job, 
                flowfile_to_read,
            )
        )

    def remove(self, pkl_too=False):
        for step in self.list_of_steps:
            step.remove()

        if pkl_too:
            os.system('rm -rf *.pkl')
                     
#endregion