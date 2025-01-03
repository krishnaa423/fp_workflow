#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
from fp.inputs.wfngeneral import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EpsilonInput:
    def __init__(
        self,
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict
        
    def get_qgrid_str(self, wfn_input: WfnGeneralInput, qshift):
        return wfn_input.get_kgrid_eps(qshift=qshift)
#endregion
