#region: Modules.
import numpy as np 
from fp.schedulers import *
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
        kpath, 
        nbands,
        job_desc,
        job_pw2bgw_desc,
    ):
        self.kpath: KPath = kpath
        self.nbands: int = nbands 
        self.job_desc: JobProcDesc = job_desc
        self.job_pw2bgw_desc: JobProcDesc = job_pw2bgw_desc
        
    def get_kgrid_str(self):
        # output = ''
        # output += 'K_POINTS crystal\n'
        # output += f'{self.kpath.shape[0]}\n'
        
        # for row in self.kpath:
        #     output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} 1.0\n'

        special_points = get_special_points(self.kpath.atoms.cell)

        output = ''
        output += 'K_POINTS {crystal_b}\n'
        output += f'{len(self.kpath.path_special_points)}\n'

        for path_special_point in self.kpath.path_special_points:
            coord = special_points[path_special_point]
            output += f'{coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f} {self.kpath.path_segment_npoints} !{path_special_point}\n'
        
        return output 

#endregion
