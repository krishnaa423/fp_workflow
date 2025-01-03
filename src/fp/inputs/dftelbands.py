#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
from fp.structure.kpath import KPath
from ase.dft.kpoints import get_special_points
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DftelbandsInput:
    def __init__(
        self,
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict 
        
    def get_kgrid(self, kpath: KPath):
        # TODO: In case of supercell and supercell unfolding. 
        output = []
        special_points = get_special_points(kpath.atoms.cell)

        for path_special_point in kpath.path_special_points:
            coord = special_points[path_special_point]
            output.append([
                float(coord[0]),
                float(coord[1]),
                float(coord[2]),
                kpath.path_segment_npoints,
                path_special_point
            ])
        
        return output 
#endregion
