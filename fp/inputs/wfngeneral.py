#region: Modules.
import numpy as np 
from fp.schedulers import *
import os
from fp.inputs.atoms import AtomsInput
from fp.structure.kpts import Kgrid 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class WfnGeneralInput:
    def __init__(
        self,
        atoms,
        kdim,
        qshift,
        is_reduced,     # Should be false for now, since using full grid.
        bands, 
        job_wfn_desc,
        job_pw2bgw_desc=None,
        job_parabands_desc=None,
        parabands_bands=None,
    ):      
        self.atoms: AtomsInput = atoms 
        self.kdim = kdim
        self.qshift = qshift 
        self.is_reduced: bool = is_reduced
        self.bands: int = bands  
        self.job_wfn_desc: JobProcDesc = job_wfn_desc
        self.job_pw2bgw_desc: JobProcDesc = job_pw2bgw_desc
        self.job_parabands_desc: JobProcDesc = job_parabands_desc
        self.parabands_bands = parabands_bands 
        
        self.kpts: np.ndarray = None 
        self.create_kgrid() 
      
    def create_kgrid(self):
        kgrid = Kgrid(
            self.atoms, 
            kdim=self.kdim,
            qshift=self.qshift,
            is_reduced=self.is_reduced
        )

        self.kpts = kgrid.get_kpts()
                    
    def get_kgrid_dft_string(self):
        output = ''
        output += 'K_POINTS crystal\n'
        
        num_kpts = self.kpts.shape[0]
        output += f'{num_kpts}\n'
        
        for row in self.kpts:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} {row[3]:15.10f}\n'
        
        return output 
   
    def get_kgrid_eps_string(self, qshift=None):
        if not qshift: qshift = self.qshift
        
        output = ''
        output += 'begin qpoints\n'
        
        # qshift.
        output += f'{qshift[0]:15.10f} {qshift[1]:15.10f} {qshift[2]:15.10f} 1.0 1\n' 
        
        # qgrid. 
        for row_idx, row in enumerate(self.kpts):
            if row_idx==0: continue 
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} 1.0 0\n'
        
        output += 'end\n'
        
        return output 
    
    def get_kgrid_sig_string(self):
        output = ''
        output += 'begin kpoints\n'
        
        # kgrid. 
        for row in self.kpts:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} 1.0\n'
        
        output += 'end\n'
        
        return output 
#endregion
