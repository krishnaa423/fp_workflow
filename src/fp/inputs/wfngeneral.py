#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
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
        input_dict: dict,
        atoms_input: AtomsInput=None,
        wfn_type: str = 'wfn',
    ):
        self.input_dict: dict = input_dict     
        
        # Extract attributes used across class methods.
        self.kdim: np.ndarray = np.array(self.input_dict[wfn_type]['kdim'])
        self.is_reduced: bool = self.input_dict[wfn_type]['sym']
        self.atoms_input: AtomsInput = atoms_input
        if atoms_input is None:
            self.atoms_input = AtomsInput(
                input_dict=self.input_dict,
            )
        if 'qshift' in self.input_dict[wfn_type]:
            self.qshift = np.array(self.input_dict[wfn_type]['qshift'])
        else:
            self.qshift = np.array([0.0, 0.0, 0.0])
        
        # Get kpoints. 
        self.kpts: np.ndarray = None 
        self.create_kgrid() 
      
    def create_kgrid(self):
        kgrid = Kgrid(
            self.atoms_input, 
            kdim=self.kdim,
            qshift=self.qshift,
            is_reduced=self.is_reduced
        )

        self.kpts = kgrid.get_kpts()
                    
    def get_kgrid_dft(self):
        output = []

        for row in self.kpts:
            output.append([
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
            ])

        return output
   
    def get_kgrid_eps(self, qshift=None):
        if not qshift: qshift = self.qshift
        
        output = []
        
        for row_idx, row in enumerate(self.kpts):
            if row_idx==0:
                output.append([
                    float(qshift[0]),
                    float(qshift[1]),
                    float(qshift[2]),
                    float(1.0),
                    int(1),
                ])
            else:
                output.append([
                    float(row[0]),
                    float(row[1]),
                    float(row[2]),
                    float(1.0),
                    int(0),
                ])
        
        return output 
    
    def get_kgrid_sig(self):
        output = []

        for row in self.kpts:
            output.append([
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(1.0),
            ])

        return output
#endregion
