#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PhmodesInput:
    def __init__(
        self,
        qidx,
        job_desc,
        extra_args: str=None,
    ):
        '''
        qidx: int 
            Starts from 1. It is the index of the irreducibe q-point.  
        '''
        self.qidx: int = qidx 
        self.job_desc: JobProcDesc = job_desc 
        self.extra_args: str = extra_args
#endregion
