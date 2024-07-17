#region: Modules.
import numpy as np 
from fp.schedulers import *
from fp.inputs.wfngeneral import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class SigmaInput:
    def __init__(
        self,
        bands,
        band_min,
        band_max,
        cutoff,
        wfn_inner_link,
        job_desc,
    ):
        self.bands = bands 
        self.band_min = band_min  
        self.band_max = band_max 
        self.cutoff = cutoff
        self.wfn_inner_link: str = wfn_inner_link
        self.job_desc: JobProcDesc = job_desc 
        
    def get_kgrid_str(self, wfn_input: WfnGeneralInput):
        return wfn_input.get_kgrid_sig_string()
#endregion
