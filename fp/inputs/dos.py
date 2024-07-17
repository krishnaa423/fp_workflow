#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DosInput:
    def __init__(
        self,
        kdim,
        bands,
        job_desc
    ):
        self.kdim: np.ndarray = np.array(kdim)
        self.bands: int = bands
        self.job_desc: JobProcDesc = job_desc
#endregion
