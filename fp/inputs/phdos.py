#region: Modules.
import numpy as np 
from fp.schedulers import *
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
    ):
        self.qdim: np.ndarray = np.array(qdim)
        self.job_desc: JobProcDesc = job_desc
#endregion
