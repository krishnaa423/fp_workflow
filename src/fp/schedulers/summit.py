#region: Modules.
from fp.schedulers.scheduler import Scheduler, JobProcDesc
from fp.io.strings import write_str_2_f
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Summit(Scheduler):
    def __init__(
        self,
        is_interactive: bool = False,
        mpi_exec: str = None,
        queue: str = None,
    ):
        super().__init__(is_interactive=is_interactive, mpi_exec=mpi_exec, queue=queue)

    def get_sched_header(self, job_desc: JobProcDesc):
        queue_str = 'batch' if not self.queue else self.queue

        # Time here in XX:xx format, where XX is the hours and xx is the minutes. 
        output = f'''#BSUB -P cph156
#BSUB -q {queue_str}
#BSUB -J struct_job
#BSUB -nnodes {job_desc.nodes}
#BSUB -W {job_desc.time}
'''
        return '\n' if self.is_interactive else output

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        return f'{self.mpi_exec} -n{job_desc.ntasks} -a1 -c7 -g1 -bpacked:7 -EOMP_NUM_THEADS=28 ' if self.mpi_exec else f'jsrun -n{job_desc.ntasks} -a1 -c7 -g1 -bpacked:7 -EOMP_NUM_THEADS=28 '
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc, add_nk_if_present: bool=True, add_ni_if_present: bool=True):
        ni = '' if not job_desc.ni or not add_ni_if_present else f' -ni {job_desc.ni} '
        nk = '' if not job_desc.nk or not add_nk_if_present else f' -nk {job_desc.nk} '
        
        output = f' {ni} {nk} '
        
        return output 

    def get_sched_submit(self, job_desc: JobProcDesc):
        return '' if self.is_interactive else 'bsub ' 
    
    def create_interactive(self, job_desc: JobProcDesc):
        string = f'''#!/bin/bash
bsub -P cph156 -q batch -J struct_job -Is -nnodes {job_desc.nodes} -W {job_desc.time} /bin/bash 
'''
        write_str_2_f('job_interactive.sh', string)
        
#endregion
