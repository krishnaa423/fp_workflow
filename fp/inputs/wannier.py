#region: Modules.
import numpy as np 
from fp.schedulers import *
from fp.inputs.atoms import *
from ase.dft.kpoints import get_special_points
from ase.data import chemical_symbols
import os 
from fp.structure import Kgrid
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class WannierInput:
    def __init__(
        self,
        atoms,
        kdim,
        num_bands,
        num_wann,
        job_wfnwan_desc,
        job_pw2wan_desc,
        job_wan_desc,
    ):
        self.atoms: AtomsInput = atoms 
        self.kdim: np.ndarray = np.array(kdim)
        self.num_bands: int = num_bands
        self.num_wann: int = num_wann
        self.job_wfnwan_desc: JobProcDesc = job_wfnwan_desc
        self.job_pw2wan_desc: JobProcDesc = job_pw2wan_desc
        self.job_wan_desc: JobProcDesc = job_wan_desc
        
    def get_unit_cell_cart(self):
        output = ''
        output += 'begin unit_cell_cart\nAng\n'
        
        for row in self.atoms.atoms.get_cell():
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
            
        output += 'end unit_cell_cart\n'
        
        return output 
    
    def get_atoms_cart(self):
        numbers = self.atoms.atoms.get_atomic_numbers()
        positions = self.atoms.atoms.get_positions()
        
        output = ''
        output += 'begin atoms_cart\nAng\n'
        
        for number, pos in zip(numbers, positions):
            output += f'{chemical_symbols[number]} {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}\n'
        
        output += 'end atoms_cart\n'
        
        return output 
        
    def get_mpgrid(self):
        output = f'mp_grid = {int(self.kdim[0])} {int(self.kdim[1])} {int(self.kdim[2])}\n'
        
        return output 
    
    def get_kpoints(self):
        kgrid = Kgrid(
            atoms=self.atoms,
            kdim=self.kdim,
            is_reduced=False,
        )

        kpts = kgrid.get_kpts()

        output = ''
        output += 'begin kpoints\n'

        for row in kpts:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} {row[3]:15.10f}\n'
        
        output += 'end kpoints\n'
        
        return output 
    
    def get_kpoints_qe(self):
        kgrid = Kgrid(
            atoms=self.atoms,
            kdim=self.kdim,
            is_reduced=False,
        )

        kpts = kgrid.get_kpts()

        output = 'K_POINTS crystal\n'
        output += f'{kpts.shape[0]}\n'

        for row in kpts:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} {row[3]:15.10f} 1.0\n'
        
        return output 
#endregion
