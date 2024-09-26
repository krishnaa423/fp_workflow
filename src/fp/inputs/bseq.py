#region: Modules.
import numpy as np 
from fp.schedulers import *
from fp.inputs.atoms import *
from fp.inputs.kernel import *
from fp.inputs.abs import *
from fp.structure.kpts import Kgrid
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class BseqInput:
    def __init__(
        self,
        atoms,
        Qdim,
        job_desc
    ):
        self.atoms: AtomsInput = atoms
        self.Qdim: np.ndarray = np.array(Qdim).astype(dtype='i4')
        self.job_desc: JobProcDesc = job_desc

    def get_Qpts(self):
        return Kgrid(self.atoms, self.Qdim, is_reduced=False).get_kpts()
#endregion
