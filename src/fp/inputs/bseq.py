#region: Modules.
import numpy as np 
from fp.inputs.atoms import AtomsInput
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
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict

    def get_Qpts(self, atoms_input: AtomsInput):
        return Kgrid(atoms_input, self.input_dict['bseq']['Qdim'], is_reduced=False).get_kpts()
#endregion
