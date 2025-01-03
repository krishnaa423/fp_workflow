#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class ScfInput:
    def __init__(
        self,
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict
        self.kdim: np.ndarray = np.array(self.input_dict['scf']['kdim'], dtype='i4')

    def get_kgrid(self):
        output = ''
        output += 'K_POINTS automatic\n'
        output += f'{int(self.kdim[0])} {int(self.kdim[1])} {int(self.kdim[2])} 0 0 0\n'
        
        return output 
#endregion
