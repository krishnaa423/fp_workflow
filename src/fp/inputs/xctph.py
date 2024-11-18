#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPhInput:
    def __init__(
        self,
        epw_qgrid,
        num_epw_val_bands,
        num_epw_cond_bands,
        num_exciton_states,
        job_desc,
    ):
        self.job_desc: JobProcDesc = job_desc
        self.num_epw_qpts: list = int(np.prod(np.array(epw_qgrid).astype('i4')))
        self.num_epw_val_bands: int = num_epw_val_bands
        self.num_epw_cond_bands: int = num_epw_cond_bands
        self.num_exciton_states: int = num_exciton_states
#endregion
