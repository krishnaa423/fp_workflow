#region: Modules.
from ase import Atoms 
from ase.optimize import BFGS 
from ase.calculators.calculator import Calculator, all_changes
from ase.data import atomic_numbers
from ase.io.trajectory import Trajectory
from ase.io.xsf import write_xsf
from ase.io import read, write
import numpy as np 
import os 
import h5py 
import xml.etree.ElementTree as ET 
import re 
import subprocess 
#endregion

#region: Variables.
iter = 0 
#endregion

#region: Functions.
def main():
    atoms = read('./sc_atoms.xsf')
    atoms.calc = CustomCalc()
    
    relax = BFGS(atoms, logfile='esd_log.txt', restart='esd_restart.pkl', trajectory='esd.traj', maxstep=15)
    relax.run(fmax=0.05)
    
    os.system('rm -rf esd_atoms.xsf')
#endregion

#region: Classes.
class CustomCalc(Calculator):
    implemented_properties = ['energy', 'forces']
    
    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):
        global iter 
        
        Calculator.calculate(self, atoms, properties, system_changes)
        
        energy, forces = self.get_energy_forces(atoms)
        
        self.results = {
            'energy': energy,
            'forces': forces,
        }
        
        iter += 1 
        
    def get_dft_energy(self):
        
        root = ET.parse('./scf.xml').getroot()
    
        elements = root.findall('.//total_energy/etot')
        
        return float(elements[0].text)*27.2114      # Energy in eV. 
    
    def get_bse_energy(self):
        '''
        TODO: Can change the exciton state being minimized. 
        '''
        return np.loadtxt('./eigenvalues.dat', skiprows=4)[0, 0]
        
    def get_energy_forces(self, atoms):
        global iter 
        
        # Create files necessary.
        write_xsf('esd_atoms.xsf', images=[atoms])
        os.system('python3 workflow.py --create')
        
        # Create savestep update. 
        with open('./job_savestep.sh', 'r') as f: txt = f.read()
        txt = re.sub(r'folder_name=.*\n', f'folder_name="./esd_save/iter_{iter}"\n', txt)
        with open('./job_savestep.sh', 'w') as f: f.write(txt)
        
        if not os.path.exists('./esd_save'): os.mkdir('./esd_save')
        if not os.path.exists(f'./esd_save/iter_{iter}'): os.mkdir(f'./esd_save/iter_{iter}')
        
        # Run once. 
        ps_result = subprocess.run(['./job_all.sh', '&>', 'job_all_esd.out'])
        
        if ps_result.returncode == 0:
            print(f'Done with iteration {iter} \n\n')
        else:
            print(f'Error in iteration {iter} with code: {ps_result.returncode}')
            os._exit(ps_result.returncode)
        
        # Get the energy and forces data. 
        energy = self.get_dft_energy() + self.get_bse_energy()
        
        with h5py.File('esf.h5', 'r') as r:
            forces = r['excited_force'][:]   
        
        return energy, forces   # eV, eV/A. 
        
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion