#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EsdStepInput:
    def __init__(
        self,
        fmax=0.1,
        max_steps=10,
        job_desc=JobProcDesc(),
    ):
        self.fmax = fmax
        self.max_steps = max_steps
        self.job_desc: JobProcDesc = job_desc
#endregion
