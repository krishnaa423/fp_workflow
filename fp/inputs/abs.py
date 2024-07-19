#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class AbsorptionInput:
    def __init__(
        self,
        val_bands_coarse,
        cond_bands_coarse,
        val_bands_fine,
        cond_bands_fine,
        Qshift,
        wfn_co_link,
        wfnq_co_link,
        wfn_fi_link,
        wfnq_fi_link,
        num_evec,
        pol_dir,
        job_desc,
    ):
        self.val_bands_coarse = val_bands_coarse  
        self.cond_bands_coarse = cond_bands_coarse 
        self.val_bands_fine = val_bands_fine 
        self.cond_bands_fine = cond_bands_fine 
        self.Qshift: np.ndarray = np.array(Qshift)
        self.wfn_co_link: str = wfn_co_link
        self.wfnq_co_link: str = wfn_co_link
        self.wfn_fi_link: str = wfn_fi_link
        self.wfnq_fi_link: str = wfnq_fi_link
        self.num_evec: int = num_evec
        self.pol_dir: np.ndarray = np.array(pol_dir)
        self.job_desc: JobProcDesc = job_desc

class PlotxctInput:
    def __init__(
        self,
        hole_position,
        supercell_size,
        state,
        wfn_fi_link,
        wfnq_fi_link,
        job_desc,
    ):
        self.hole_position = hole_position 
        self.supercell_size = supercell_size
        self.state = state 
        self.wfn_fi_link: str = wfn_fi_link
        self.wfnq_fi_link: str = wfnq_fi_link
        self.job_desc: JobProcDesc = job_desc
        
    def get_hole_position_str(self):
        output = f'{self.hole_position[0]:15.10f} {self.hole_position[1]:15.10f} {self.hole_position[2]:15.10f}'
        
        return output 
    
    def get_supercell_size_str(self):
        output = f'{int(self.supercell_size[0])} {int(self.supercell_size[1])} {int(self.supercell_size[2])}'
        
        return output 
#endregion
