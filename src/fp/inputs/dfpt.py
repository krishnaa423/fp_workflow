#region: Modules.
from fp.inputs.atoms import *
from fp.schedulers.scheduler import *
import numpy as np
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DfptInput:
    def __init__(
        self,
        atoms,
        qgrid,
        conv_threshold,
        job_desc,
        extra_args: str=None,
    ):
        self.atoms: AtomsInput = atoms
        self.qgrid: np.array = np.array(qgrid)
        self.conv_threshold: float = conv_threshold
        self.job_desc: JobProcDesc = job_desc
        self.extra_args: str = extra_args
#endregion