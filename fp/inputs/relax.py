#region: Modules.
import numpy as np
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class RelaxType:
    GS_RELAX = 0
    GS_VC_RELAX = 1
    CDFT_RELAX = 2
    CDFT_VC_RELAX = 3

class RelaxInput:
    def __init__(
        self, 
        max_val: int,
        job_desc,
        relax_type: int = RelaxType.GS_RELAX,
    ):
        self.job_desc: JobProcDesc = job_desc
        self.max_val: int = max_val
        self.relax_type: RelaxType = relax_type
        
    def get_occupations_str(self):

        nbands = self.max_val if self.relax_type==RelaxType.GS_RELAX or self.relax_type==RelaxType.GS_VC_RELAX else self.max_val+1

        occupations = np.zeros((nbands,), dtype='f8')

        if self.relax_type==RelaxType.GS_RELAX or self.relax_type==RelaxType.GS_VC_RELAX:
            occupations[:] = 2.0
        elif self.relax_type==RelaxType.CDFT_RELAX or self.relax_type==RelaxType.CDFT_VC_RELAX:
            occupations[:-2] = 2.0
            occupations[-1] = 2.0
        else:
            Exception('Must set valid occupations type')

        output = 'OCCUPATIONS \n'
        
        for occ in occupations:
            output += f'{occ:15.10f}\n'
            
        return output 
    
    def get_nbnd(self):
        nbands = self.max_val if self.relax_type==RelaxType.GS_RELAX or self.relax_type==RelaxType.GS_VC_RELAX else self.max_val+1

        return int(nbands)
    
    def calc_str(self):
        output = 'relax' if self.relax_type==RelaxType.GS_RELAX or self.relax_type==RelaxType.CDFT_RELAX else 'vc-relax'

        return output
#endregion