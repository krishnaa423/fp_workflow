#region: Modules.
import numpy as np 
from fp.schedulers import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EpwInput:
    def __init__(
        self,
        kgrid_coarse,
        qgrid_coarse,
        kgrid_fine,
        qgrid_fine,
        bands,
        job_desc,
        exec_loc,
        skipped_bands=None,
    ): 
        self.kgrid_coarse:np.ndarray = kgrid_coarse
        self.qgrid_coarse:np.ndarray = qgrid_coarse
        self.kgrid_fine:np.ndarray = kgrid_fine
        self.qgrid_fine:np.ndarray = qgrid_fine
        self.bands = bands  
        self.job_desc: JobProcDesc = job_desc
        self.exec_loc: str = exec_loc
        self.skipped_bands: list[tuple] = skipped_bands
        
    def get_skipped_bands_str(self):
        bands_skipped = self.skipped_bands
        
        bands_skipped_str = ''
        
        if bands_skipped:
            num_bands_skipped = len(bands_skipped)
            exclude_bands_str = "'exclude_bands="
            
            for bands_idx, bands in enumerate(bands_skipped):
                exclude_bands_str += f'{bands[0]}:{bands[1]}'
                if bands_idx!=num_bands_skipped-1: exclude_bands_str += ','
                
            exclude_bands_str += "'"
            
            bands_skipped_str = 'bands_skipped=' + exclude_bands_str
        
        return bands_skipped_str
#endregion
