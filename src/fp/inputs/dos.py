#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DosInput:
    def __init__(
        self,
        kdim,
        bands,
        job_desc,
        extra_control_args: str=None,
        extra_system_args: str=None,
        extra_electron_args: str=None,
        extra_dos_args: str=None,
        extra_pdos_args: str=None,
    ):
        self.kdim: np.ndarray = np.array(kdim)
        self.bands: int = bands
        self.job_desc: JobProcDesc = job_desc
        self.extra_control_args: str = extra_control_args
        self.extra_system_args: str = extra_system_args
        self.extra_electron_args: str = extra_electron_args
        self.extra_dos_args: str = extra_dos_args
        self.extra_pdos_args: str = extra_pdos_args
#endregion
