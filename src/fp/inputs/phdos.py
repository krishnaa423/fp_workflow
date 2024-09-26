#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PhdosInput:
    def __init__(
        self,
        qdim,
        job_desc,
        extra_q2r_args: str=None,
        extra_matdyn_args: str=None,
    ):
        self.qdim: np.ndarray = np.array(qdim)
        self.job_desc: JobProcDesc = job_desc
        self.extra_q2r_args: str = extra_q2r_args
        self.extra_matdyn_args: str = extra_matdyn_args
#endregion
