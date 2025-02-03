#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EpwInput:
    def __init__(
        self,
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict
        
    def get_skipped_bands_str(self, bands_skipped: list=None):
        # Populate list.
        if bands_skipped is None:
            bands_skipped = []
            abs_val_bands = self.input_dict['abs']['num_val_bands']
            total_val_bands = self.input_dict['total_valence_bands']
            abs_cond_bands = self.input_dict['abs']['num_cond_bands']
            wfn_cond = self.input_dict['wfn']['num_cond_bands']

            if abs_val_bands!= total_val_bands:
                temp = (1, total_val_bands - abs_val_bands)
                bands_skipped.append(temp)

            if abs_cond_bands!= wfn_cond and abs_cond_bands < wfn_cond:
                temp = (total_val_bands+ abs_cond_bands + 1, wfn_cond + total_val_bands)
                bands_skipped.append(temp)

            if len(bands_skipped)==0:
                bands_skipped = None

        # Populate string. 
        bands_skipped_str = ''
        exclude_bands_str = None
        if bands_skipped is not None:
            num_bands_skipped = len(bands_skipped)
            exclude_bands_str = "'exclude_bands="
            
            for bands_idx, bands in enumerate(bands_skipped):
                exclude_bands_str += f'{bands[0]}:{bands[1]}'
                if bands_idx!=num_bands_skipped-1: exclude_bands_str += ','
                
            exclude_bands_str += "'"
            
            bands_skipped_str = 'bands_skipped=' + exclude_bands_str
        
        return exclude_bands_str
#endregion
