#region: Modules.
from fp.schedulers.scheduler import Scheduler, JobProcDesc
from fp.io.strings import write_str_2_f
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Frontera(Scheduler):
    def __init__(
        self,
        is_interactive: bool = False,
        queue: str = None,
    ):
        '''
        Queue could be development or normal.
        '''
        super().__init__(is_interactive)
        self.queue: str = queue

    def get_sched_header(self, job_desc: JobProcDesc):
        partition_str = 'development' if not self.queue else self.queue

        # Time here in XX:xx:yy format, where XX is the hours and xx is the minutes, yy is seconds.
        output = f'''#SBATCH --account=PHY20032
#SBATCH --partition={partition_str}
#SBATCH --job-name=struct_job
#SBATCH --nodes={job_desc.nodes}
#SBATCH --time={job_desc.time}
#SBATCH --mail-type=all
$SBATCH --mail-user=krishnaa.vadivel@yale.edu
'''
        return '\n' if self.is_interactive else output

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        output = f'ibrun --ntasks={job_desc.ntasks} ' 
        
        return output 
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc):
        ni = '' if not job_desc.ni else f' -ni {job_desc.ni} '
        nk = '' if not job_desc.nk else f' -nk {job_desc.nk} '
        
        output = f' {ni} {nk} '
        
        return output 

    def get_sched_submit(self):
        return '' if self.is_interactive else 'sbatch ' 
    
    def create_interactive(self, job_desc: JobProcDesc):
        string = f'''#!/bin/bash
srun --account=PHY20032 --partition=normal --pty --job-name=struct_job --nodes={job_desc.nodes} --time={job_desc.time} /bin/bash -l
'''
        write_str_2_f('job_interactive.sh', string)
        
#endregion
