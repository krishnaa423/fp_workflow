#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DftelbandsInput:
    def __init__(
        self,
        kpath, 
        nbands,
        job_desc,
        job_pw2bgw_desc,
    ):
        self.kpath:np.ndarray = np.array(kpath) 
        self.nbands: int = nbands 
        self.job_desc: JobProcDesc = job_desc
        self.job_pw2bgw_desc: JobProcDesc = job_pw2bgw_desc
        
    def get_kgrid_str(self):
        output = ''
        output += 'K_POINTS crystal\n'
        output += f'{self.kpath.shape[0]}\n'
        
        for row in self.kpath:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} 1.0\n'
        
        return output 

#endregion
