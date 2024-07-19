#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class KernelInput:
    def __init__(
        self,
        val_bands_coarse,
        cond_bands_coarse,
        Qshift,
        wfn_co_link,
        wfnq_co_link,
        job_desc,
    ):
        self.val_bands_coarse = val_bands_coarse  
        self.cond_bands_coarse = cond_bands_coarse 
        self.Qshift: np.ndarray = np.array(Qshift)
        self.wfn_co_link: str = wfn_co_link
        self.wfnq_co_link: str = wfnq_co_link
        self.job_desc: JobProcDesc = job_desc
#endregion
