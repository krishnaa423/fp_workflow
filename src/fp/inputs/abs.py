#region: Modules.
import numpy as np 
from fp.schedulers.scheduler import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PlotxctInput:
    def __init__(self, input_dict:dict):
        self.input_dict: dict = input_dict
        
    def get_hole_position_str(self):
        output = f'{self.input_dict['plotxct']['hole_position'][0]:15.10f} {self.input_dict['plotxct']['hole_position'][1]:15.10f} {self.input_dict['plotxct']['hole_position'][2]:15.10f}'
        
        return output 
    
    def get_supercell_size_str(self):
        output = f'{int(self.input_dict['plotxct']['supercell_size'][0])} {int(self.input_dict['plotxct']['supercell_size'][1])} {int(self.input_dict['plotxct']['supercell_size'][2])}'
        
        return output 
#endregion
