#region: Modules.
from ase.units import Hartree
from ase.calculators.calculator import Calculator, all_changes
import numpy as np 
from fp.analysis.dft import *
from fp.analysis.elph import *
from fp.analysis.eqp import *
from fp.analysis.xct import *
from fp.analysis.esf import *
from ase.io.xsf import write_xsf
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EsfCalculator(Calculator):
    implemented_properties = ['energy', 'forces']
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.excited_force = None 
        
    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):
        super().calculate(atoms, properties=properties, system_changes=system_changes)
        
        # TODO. 
        self.results['energy'] = 0.0

        with h5py.File('esf.h5', 'r') as r:
            self.excited_force: np.ndarray = r['excited_force'][:]

        # self.results['forces'] = self.excited_force.flatten()*Hartree  # I don't think multiplying by Hartree is correct. 
        self.results['forces'] = self.excited_force

class Esf:
    def __init__(self):
        '''
        Class object variables declared and defined in self.assemble_components() function. 
        '''
        pass 
    
    def assemble_components(self):
        self.atoms = DftAtoms().get_atoms()
        
        self.dft_force = DftForceResult().get_dftforce()
        
        self.xct = XctResult()
        self.vbm, self.nc, self.nv = self.xct.get_vbm_nc_nv()
        
        self.eig_c, self.eig_v = EqpResult().get_eig(self.vbm, self.nc, self.nv)
        
        self.elph_c, self.elph_v = ElphResult().get_elph(self.vbm, self.nc, self.nv)
        
        self.xct_evec = self.xct.get_xctevec()
    
    def do_calculation(self):
        F = np.zeros(shape=(self.atoms.get_number_of_atoms()*3,), dtype='c16')
        
        
        # Sum as per method in Strubbe's thesis. page 162. 
        F += np.einsum(
            'scc,cv,cv->s',
            self.elph_c,
            np.conjugate(self.xct_evec),
            self.xct_evec,
        )
        
        F -= np.einsum(
            'svv,cv,cv->s',
            self.elph_v,
            np.conjugate(self.xct_evec),
            self.xct_evec,
        )
        
        F += np.einsum(
            'Cv,scC,cC,cv->s',
            np.conjugate(self.xct_evec),
            self.elph_c,
            self.eig_c,
            self.xct_evec,
        )
        
        F += np.einsum(
            'CV,sVv,Vv,cv->s',
            np.conjugate(self.xct_evec),
            self.elph_v,
            self.eig_v,
            self.xct_evec,
        )
        
        self.excited_force = np.real(F).reshape(int(F.size/3), 3) + self.dft_force
           
    def calc(self):
        self.assemble_components()
        self.do_calculation()
    
    def write(self):
        
        # Write hdf5 file. 
        with h5py.File('esf.h5', 'w') as f:
            natoms = self.atoms.get_number_of_atoms()
            ds_positions = f.create_dataset('positions', shape=(natoms, 3))
            ds_dft_force = f.create_dataset('dft_force', shape=(natoms, 3))
            ds_excited_force = f.create_dataset('excited_force', shape=(natoms, 3))
            
            ds_positions[:] = self.atoms.get_positions()
            ds_dft_force[:] = self.dft_force
            ds_excited_force[:] = self.excited_force
            
        # Write xsf file. 
        self.atoms.calc = EsfCalculator()
        self.atoms.get_potential_energy()
        self.atoms.get_forces()
        with open('esf.xsf', 'w') as w: write_xsf(w, [self.atoms])
#endregion
