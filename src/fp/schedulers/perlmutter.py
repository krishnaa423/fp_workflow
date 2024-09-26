#region: Modules.
from fp.schedulers.scheduler import Scheduler, JobProcDesc
from fp.io.strings import write_str_2_f
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Perlmutter(Scheduler):
    def __init__(
        self,
        is_interactive: bool = False,
        mpi_exec: str = None,
        queue: str = None,
        is_gpu: bool = False,
    ):
        super().__init__(is_interactive=is_interactive, mpi_exec=mpi_exec, queue=queue)
        self.is_gpu: bool = is_gpu

        self.constraint_str = 'gpu' if self.is_gpu else 'cpu'

    def get_sched_header(self, job_desc: JobProcDesc):
        qos_str = 'debug' if not self.queue else self.queue

        # Time here in XX:xx format, where XX is the hours and xx is the minutes. 
        output = f'''#SBATCH --account=m3571
#SBATCH --qos={qos_str}
#SBATCH --job-name=struct_job
#SBATCH --constraint={self.constraint_str}
#SBATCH --nodes={job_desc.nodes}
#SBATCH --time={job_desc.time}
#SBATCH --mail-type=all
$SBATCH --mail-user=krishnaa.vadivel@yale.edu
'''
        return '\n' if self.is_interactive else output

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        output = f'{self.mpi_exec} --ntasks={job_desc.ntasks} ' if self.mpi_exec else f'srun --ntasks={job_desc.ntasks} '
        if self.is_gpu: output += '  --gpus-per-task=1 '  # --cpus-per-task=16 with OMP_NUM_THREADS=2.
        
        return output 
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc, add_nk_if_present: bool=True, add_ni_if_present: bool=True):
        ni = '' if not job_desc.ni or not add_ni_if_present else f' -ni {job_desc.ni} '
        nk = '' if not job_desc.nk or not add_nk_if_present else f' -nk {job_desc.nk} '
        
        output = f' {ni} {nk} '
        
        return output 

    def get_sched_submit(self):
        return '' if self.is_interactive else 'sbatch ' 
    
    def create_interactive(self, job_desc: JobProcDesc):
        string = f'''#!/bin/bash
salloc --account=m3571 --qos=interactive --job-name=struct_job --constraint={self.constraint_str} --nodes={job_desc.nodes} --time={job_desc.time}
'''
        write_str_2_f('job_interactive.sh', string)
        
#endregion
