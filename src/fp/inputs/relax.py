#region: Modules.
import numpy as np
from fp.schedulers.scheduler import *
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
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict

        # Extract some variables to use across class methods. 
        self.relax_type: str = self.input_dict['relax']['type']
        self.is_spinorbit: bool = self.input_dict['scf']['is_spinorbit']

    def get_occupations(self):
        self.total_valence_bands: int = self.input_dict['total_valence_bands']

        nbands = self.total_valence_bands+1 if 'cdft' in self.relax_type else self.total_valence_bands

        occupations = np.zeros((nbands,), dtype='f8')

        if 'cdft' in self.relax_type:
            if self.is_spinorbit:
                occupations[:-2] = 1.0
                occupations[-1] = 1.0
            else:
                occupations[:-2] = 2.0
                occupations[-1] = 2.0
        else: 
            if self.is_spinorbit:
                occupations[:] = 1.0
            else:
                occupations[:] = 2.0
            
        return occupations.tolist()
    
    def get_nbnd(self):
        self.total_valence_bands: int = self.input_dict['total_valence_bands']
        nbands = self.total_valence_bands+1 if 'cdft' in self.relax_type else self.total_valence_bands

        return int(nbands)
    
    def calc_str(self):
        output = 'relax' if self.relax_type in ['relax', 'cdft-relax'] else 'vc-relax'

        return output
#endregion