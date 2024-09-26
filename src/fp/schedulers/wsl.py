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
        mpi_exec: str = None,
        is_gpu: bool = False,
    ):
        super().__init__(is_interactive=is_interactive, mpi_exec=mpi_exec)
        self.is_gpu: bool = is_gpu

    def get_sched_header(self, job_desc: JobProcDesc):
        return ''

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        return '' if self.is_gpu else f'{self.mpi_exec} -n {job_desc.ntasks} ' if self.mpi_exec else f'mpirun -n {job_desc.ntasks} ' 
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc, add_nk_if_present: bool=True, add_ni_if_present: bool=True):
        ni = '' if not job_desc.ni or not add_ni_if_present else f' -ni {job_desc.ni} '
        nk = '' if not job_desc.nk or not add_nk_if_present else f' -nk {job_desc.nk} '
        
        output = f' {ni} {nk} '
        
        return output 

    def get_sched_submit(self):
        return '' 
        
#endregion
