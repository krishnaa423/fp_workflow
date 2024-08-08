#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPolInput:
    def __init__(
        self,
        max_error,
        max_steps,
        job_desc
    ):
        self.max_error = max_error
        self.max_steps = max_steps
        self.job_desc: JobProcDesc = job_desc
#endregion
