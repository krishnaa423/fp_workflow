#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class KpdosInput:
    def __init__(
        self,
        job_desc,
        extra_args: str=None,
    ):
        self.job_desc: JobProcDesc = job_desc
        self.extra_args: str = extra_args
#endregion
