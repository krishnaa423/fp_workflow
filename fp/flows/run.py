#region: Modules.
import subprocess
import time 
import os 
from fp.inputs.input_main import *
#endregion

#region: Variables.
#endregion

#region: Functions.
def run_and_wait_command(script, input: Input, total_time):
    '''
    Run each script and write out some logging info. 
    '''

    start_time = time.time()
    print(f'Starting {script}.', flush=True)
    ps_result = subprocess.run(f'{input.scheduler.get_sched_submit()}{script}')
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    total_time += elapsed_time

    if ps_result.returncode == 0:  # Success.
        print(f'Done with {script} in {elapsed_time:15.10f} seconds.\n\n', flush=True)
    else:               # Fail.
        print(f'Error finishing: {script}. Exited with code {ps_result.returncode}. Time elapsed is {elapsed_time:15.10f} seconds.\n\n', flush=True)
        print(f'Total time for workflow run in {total_time:15.10f} seconds.\n', flush=True)
        os._exit(ps_result.returncode)

    return total_time
#endregion

#region: Classes.
#endregion
