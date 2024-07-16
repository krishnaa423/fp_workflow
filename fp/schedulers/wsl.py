#region: Modules.
from fp.schedulers.scheduler import Scheduler, JobProcDesc
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class WSL(Scheduler):
    def __init__(
        self,
        is_interactive: bool = False,
        is_gpu: bool = False,
    ):
        super().__init__(is_interactive)
        self.is_gpu: bool = is_gpu

    def get_sched_header(self, job_desc: JobProcDesc):
        return ''

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        return '' if self.is_gpu else f'mpirun -n {job_desc.ntasks} ' 
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc):
        ni = '' if not job_desc.ni else f' -ni {self.job_desc.ni} '
        nk = '' if not job_desc.nk else f' -nk {self.job_desc.nk} '
        
        output = f' {ni} {nk} '
        
        return output 

    def get_sched_submit(self):
        return '' 
        
#endregion
