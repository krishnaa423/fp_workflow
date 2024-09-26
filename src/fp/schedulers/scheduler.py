#region: Modules.
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class JobProcDesc:
    def __init__(
        self,
        nodes: int = None,
        ntasks: int = None,
        time: int = None,
        nk: int = None,
        ni: int = None,
    ):
        self.nodes: int = nodes
        self.ntasks: int = ntasks
        self.time: str = time
        self.nk: int = nk
        self.ni: int = ni

class Scheduler:
    def __init__(
        self,
        is_interactive:bool = False,
        mpi_exec: str = None,
        queue: str = None,
    ):
        self.is_interactive: bool = is_interactive
        self.mpi_exec: str = mpi_exec
        self.queue: str = queue

    def get_sched_header(self, job_desc: JobProcDesc):
        return ''

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        return ''
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc, add_nk_if_present: bool=True, add_ni_if_present: bool=True):
        return ''

    def get_sched_submit(self):
        return '' 
#endregion
