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
class SigmaInput:
    def __init__(
        self,
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict
        
    def get_kgrid(self, wfn_input: WfnGeneralInput):
        return wfn_input.get_kgrid_sig()
#endregion
