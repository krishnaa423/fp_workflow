#region: Modules.
from ase import Atoms 
import numpy as np 
import os
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Kgrid:
    def __init__(
        self, 
        atoms,
        kdim,
        is_reduced=False,
    ):
        self.kdim: np.ndarray = np.array(kdim).astype(dtype='i4')
        self.atoms: Atoms = atoms 
        self.is_reduced: bool = is_reduced

    def get_fbz_kpts(self):
        # Calc the kpts.
        os.system(f'kmesh.pl {int(self.kdim[0])} {int(self.kdim[1])} {int(self.kdim[2])} &> kmesh.pl.out')

        # Read the kpts.
        kpts = np.loadtxt('kmesh.pl.out', skiprows=2, dtype='f8')

        # Delete the generated files.
        os.system('rm -rf kmesh.pl.out')

        return kpts

    def get_ibz_kpts(self):
        pass 

    def get_kpts(self):
        return self.get_ibz_kpts() if self.is_reduced else self.get_fbz_kpts()
#endregion
