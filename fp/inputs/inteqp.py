#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class InteqpInput:
    def __init__(
        self,
        val_bands_coarse,
        cond_bands_coarse,
        val_bands_fine,
        cond_bands_fine,
        wfn_co_link,
        wfn_fi_link,
        job_desc,
    ):
        self.val_bands_coarse = val_bands_coarse  
        self.cond_bands_coarse = cond_bands_coarse 
        self.val_bands_fine = val_bands_fine 
        self.cond_bands_fine = cond_bands_fine 
        self.wfn_co_link: str = wfn_co_link
        self.wfn_fi_link: str = wfn_fi_link
        self.job_desc: JobProcDesc = job_desc

#endregion
