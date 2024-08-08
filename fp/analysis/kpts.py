#region: Modules.
from fp.inputs.input_main import Input
from fp.io import load_obj
from fp.structure import Kgrid
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class KptsResult:
    def __init__(
        self,
        inputpkl_filename,
    ):
        self.inputpkl_filename = inputpkl_filename
        self.input: Input = load_obj(self.inputpkl_filename)

        # After class method calls. 
        self.kpts: np.ndarray = None 

    def get_wfn_kpts(self):
        kpts = Kgrid(
            self.input.atoms,
            self.input.wfn.kdim,
            qshift=self.input.wfn.qshift,
            is_reduced=self.input.wfn.is_reduced,
        )
        self.kpts = kpts.get_kpts()
        
        return self.kpts
        
#endregion
