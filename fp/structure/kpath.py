#region: Modules.
from ase import Atoms 
from ase.dft.kpoints import BandPath
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class KPath:
    def __init__(
            self, 
            atoms: Atoms, 
            path_string: str, 
            npoints: int
        ):
        self.atoms: Atoms = atoms
        self.path_string: str = path_string
        self.npoints: int = npoints

        # generate bandpath.
        self.bandpath = atoms.cell.bandpath(path=self.path_string, npoints=self.npoints)

    def get_kpts(self):
        return self.bandpath.kpts
    
    def get_axis(self):
        return self.bandpath.get_linear_kpoint_axis()
    
    def get_sc_path(self, sc_grid: np.ndarray):
        pass 
#endregion
