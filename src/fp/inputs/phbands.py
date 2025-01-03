#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
from fp.structure.kpath import KPath
from ase import Atoms
from ase.dft.kpoints import get_special_points
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PhbandsInput:
    def __init__(
        self,
        input_dict: dict,
        atoms: Atoms,
    ):
        self.input_dict: dict = input_dict
        self.kpath: KPath = KPath(self.input_dict, atoms)
        
    def get_kpath_str(self):
        output = ''
        special_points = get_special_points(self.kpath.atoms.cell)

        output += f'{len(self.kpath.path_special_points)}\n'

        for path_special_point in self.kpath.path_special_points:
            coord = special_points[path_special_point]
            output += f'{coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f} {self.kpath.path_segment_npoints} !{path_special_point}\n'
        
        return output 
#endregion
