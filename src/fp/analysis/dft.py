#region: Modules.
from ase.io import read 
import xml.etree.ElementTree as ET
import numpy as np 
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DftAtomsResult:
    def get_atoms(self):
        atoms = read('./esd_atoms.xsf') if os.path.exists('./esd_atoms.xsf') else read('./sc_atoms.xsf')
        
        return atoms 

class DftForceResult:
    '''
    Returns result in eV/A. 
    '''
    def get_dftforce(self, inevA_units=False):
        habohr2eva = 27.2114/0.529177 
        root = ET.parse('./scf.xml').getroot()
    
        elements = root.findall('.//output/forces')
        
        # Set sizes. 
        dft_forces: np.ndarray = np.fromstring(elements[0].text, dtype='f8', sep=' ')
        if inevA_units: dft_forces *= habohr2eva     # Ha/bohr -> eV/A. 
        dft_forces = dft_forces.reshape(int(dft_forces.size/3), 3)
        
        return dft_forces
#endregion
