#region: Modules.
from fp.inputs.atoms import *
from fp.schedulers import *
import numpy as np
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DfptInput:
    def __init__(
        self,
        atoms,
        qgrid,
        conv_threshold,
        job_desc,
        njobs=None,
    ):
        self.atoms: AtomsInput = atoms
        self.qgrid: np.array = np.array(qgrid)
        self.conv_threshold: float = conv_threshold
        self.njobs: int = njobs
        self.job_desc: JobProcDesc = job_desc
    
    def get_irr_lim(self, job_idx):
        nmodes: int = self.atoms.atoms.get_number_of_atoms()*3
        
        if not self.njobs or self.njobs==1:
            return 1, nmodes
        
        if self.njobs <= nmodes:
            
            modes_per_job = nmodes // self.njobs  
            modes_per_job_last = nmodes - modes_per_job*(self.njobs  - 1)
            
            start_irr = job_idx*modes_per_job + 1 
            last_irr = start_irr + modes_per_job - 1  if job_idx != self.njobs-1 else start_irr + modes_per_job_last - 1 
            
        else:
            raise Exception('njobs <= nmodes should be satisfied.')
        
        return start_irr, last_irr 
#endregion