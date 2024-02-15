import os 
import subprocess
import json 
import time 

# region: Varibles
total_time: float = 0.0
input: dict = None 
run_order:list = [
    "./job_scf.sh",
    
    "./job_wfn.sh",
    "./job_wfn_pw2bgw.sh",
    "./job_parabands.sh",
    
    "./job_wfnq.sh",
    "./job_wfnq_pw2bgw.sh",
    
    # "./job_wfn_fi.sh",
    # "./job_abi2bgw.sh",
    
    # "./job_elband.sh",
    
    # "./job_fold2bloch.sh",
    
    # "job_wannier.sh",
    
    # "job_gwbse.sh",
    
    # "job_phband.sh",
    # "job_mrgddb.sh",
    # "job_anaddb.sh",
    
    "./job_epsilon.sh",
    
    # "./job_sigma.sh",
    
    # "./job_kernel.sh",
    
    # "./job_absorption.sh",
    # "./job_plotxct.sh",
    
    # "job_esf.sh",
    
    # "job_step.sh",
    
    # "job_update.sh",
]
# endregion

# region: Functions

def get_input():
    global input
    
    with open('input.json', 'r') as f: input = json.load(f)

def run_task(script):
    '''
    Run each script and time it. 
    Exit if there was as error. 
    '''
    global total_time 
    
    start_time = time.time()
    print(f'Starting {script}.', flush=True)
    ps_result = subprocess.run(f'{input["sched_launch"]}{script}')
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    total_time += elapsed_time
    
    if ps_result.returncode == 0:  # Success. 
        print(f'Done with {script} in {elapsed_time:15.10f} seconds.\n\n', flush=True)
    else:               # Fail. 
        print(f'Error finishing: {script}. Exited with code {ps_result.returncode}. Time elapsed is {elapsed_time:15.10f} seconds.\n\n', flush=True)
        print(f'Total time for workflow run in {total_time:15.10f} seconds.\n', flush=True)
        os._exit(ps_result.returncode)

# endregion

# region: Main
get_input()
for run in run_order:
    run_task(run)
print(f'Total time for workflow run in {total_time:15.10f} seconds.\n', flush=True)
# endregion