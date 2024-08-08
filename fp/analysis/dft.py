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
class DftAtoms:
    def get_atoms(self):
        atoms = read('./esd_atoms.xsf') if os.path.exists('./esd_atoms.xsf') else read('./sc_atoms.xsf')
        
        return atoms 

class DftForceResult:
    def get_dftforce(self):
        root = ET.parse('./scf.xml').getroot()
    
        elements = root.findall('.//output/forces')
        
        # Set sizes. 
        dft_forces = np.fromstring(elements[0].text, dtype='f8', sep=' ')*27.2114/0.529177      # Ha/bohr -> eV/A. 
        dft_forces = dft_forces.reshape(int(dft_forces.size/3), 3)
        
        return dft_forces
#endregion
