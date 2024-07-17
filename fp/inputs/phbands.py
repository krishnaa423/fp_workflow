#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PhbandsInput:
    def __init__(
        self,
        kpath,
        job_desc,
    ):
        self.kpath:np.ndarray = np.array(kpath) 
        self.job_desc: JobProcDesc = job_desc
        
    def get_kpath_str(self):
        output = ''
        output += f'{self.kpath.shape[0]}\n'
        
        for row in self.kpath:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        
        return output 
#endregion
