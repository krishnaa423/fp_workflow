#region: Modules.
from ase import Atoms 
import numpy as np 
from ase.data import chemical_symbols, atomic_masses
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class AtomsInput:
    def __init__(
        self,
        atoms,
    ):
        self.atoms: Atoms = atoms  
        
    def get_ntyp(self):
        return len(np.unique(self.atoms.get_atomic_numbers()))
    
    def get_nat(self):
        return len(self.atoms)

    def get_scf_cell(self):
        output = ''
        for row in self.atoms.get_cell():
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        return output 
    
    def get_scf_atomic_species(self):
        output = ''
        
        for atm_num in np.unique(self.atoms.get_atomic_numbers()):
            output += f'{chemical_symbols[atm_num]} {atomic_masses[atm_num]} {chemical_symbols[atm_num]}_ONCV_PBE_sr.upf\n'
        return output 

    def get_scf_atomic_positions(self, first_column='symbol'):
        output = ''
        
        if first_column=='symbol':
            for atm_num, row in zip(self.atoms.get_atomic_numbers(), self.atoms.get_positions()):
                output += f'{chemical_symbols[atm_num]} {row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        
        if first_column=='atom_index':
            _, atom_index = np.unique(self.atoms.get_atomic_numbers(), return_inverse=True)
            atom_index += 1     # 1 based index.
            for atm_num, row in zip(atom_index, self.atoms.get_positions()):
                output += f'{atm_num} {row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        return output 
#endregion
