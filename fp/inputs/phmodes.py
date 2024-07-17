#region: Modules.
import numpy as np 
from fp.schedulers import *
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
    ):
        '''
        qidx: int 
            Starts from 1. It is the index of the irreducibe q-point.  
        '''
        self.qidx: int = qidx 
        self.job_desc: JobProcDesc = job_desc 
#endregion
