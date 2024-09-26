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
class PhbandsInput:
    def __init__(
        self,
        kpath,
        job_desc,
        extra_q2r_args: str=None,
        extra_matdyn_args: str=None,
    ):
        self.kpath: KPath = kpath
        self.job_desc: JobProcDesc = job_desc
        self.extra_q2r_args: str = extra_q2r_args
        self.extra_matdyn_args: str = extra_matdyn_args
        
    def get_kpath_str(self):
        output = ''
        if self.kpath.path_total_npoints:
            kpts = self.kpath.get_kpts()
            output += f'{kpts.shape[0]}\n'
            
            for row in self.kpath:
                output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        else:
            special_points = get_special_points(self.kpath.atoms.cell)

            output += f'{len(self.kpath.path_special_points)}\n'

            for path_special_point in self.kpath.path_special_points:
                coord = special_points[path_special_point]
                output += f'{coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f} {self.kpath.path_segment_npoints} !{path_special_point}\n'
        
        return output 
#endregion
