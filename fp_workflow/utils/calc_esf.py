#region: Modules.
import h5py 
import os 
from ase import Atoms 
from ase.units import Hartree
from ase.io import read, write 
from ase.io.xsf import write_xsf
import xml.etree.ElementTree as ET 
import numpy as np 
from ase.calculators.calculator import Calculator, all_changes
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    '''
    Assume gamma point supercell for calculations. 
    '''
    esf = Esf()
    esf.assemble_components()
    esf.do_calculation()
    esf.write()
#endregion

#region: Classes.
class EsfCalculator(Calculator):
    implemented_properties = ['energy', 'forces']
    
    def __init__(self, excited_force, **kwargs):
        super().__init__(**kwargs)
        self.excited_force = excited_force*Hartree
        
    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):
        super().calculate(atoms, properties=properties, system_changes=system_changes)
        
        self.results['energy'] = 0.0
        self.results['forces'] = self.excited_force

class DftAtoms:
    def get_atoms(self):
        atoms = read('./esd_atoms.xsf') if os.path.exists('./esd_atoms.xsf') else read('./sc_atoms.xsf')
        
        return atoms 

class DftForce:
    def get_dftforce(self):
        root = ET.parse('./scf.xml').getroot()
    
        elements = root.findall('.//output/forces')
        
        # Set sizes. 
        dft_forces = np.fromstring(elements[0].text, dtype='f8', sep=' ')*27.2114/0.529177      # Ha/bohr -> eV/A. 
        dft_forces = dft_forces.reshape(int(dft_forces.size/3), 3)
        
        return dft_forces

class Eqp:
    def __init__(self):
        self.dft: np.ndarray = None  
        self.gw: np.ndarray = None 
    
    def get_dfteig(self):
        data = np.loadtxt('./bandstructure.dat', skiprows=2, dtype='f8')
        
        self.bands = data[:, 1].astype(dtype='i4')
        self.dft=  data[:, 5]
        
        return self.dft
    
    def get_gweig(self):
        data = np.loadtxt('./bandstructure.dat', skiprows=2, dtype='f8')
        
        self.bands = data[:, 1].astype(dtype='i4')
        self.gw = data[:, 6]
        
        return self.gw 
    
    def get_eig(self, vbm, nc, nv):
        self.get_dfteig()
        self.get_gweig()
        
        vbm_idx = np.where(self.bands == vbm)[0][0]
        
        dft_eig_c = self.dft[vbm_idx+1::1]
        dft_eig_v = self.dft[vbm_idx::-1]
        
        gw_eig_c = self.gw[vbm_idx+1::1]
        gw_eig_v = self.gw[vbm_idx::-1]
        
        # create cond.
        eig_c = np.zeros(shape=(nc, nc), dtype='c16')
        for i in range(nc):
            for j in range(nc):
                if i == j:
                    eig_c[i, j] = 1.0
                else:
                    den = dft_eig_c[i] - dft_eig_c[j]
                    if den < 1e-12:
                        eig_c[i, j] = 1.0
                    else:
                        num = gw_eig_c[i] - gw_eig_c[j]
                        eig_c[i, j] = num/den 
            
        # create val. 
        eig_v = np.zeros(shape=(nv, nv), dtype='c16')
        for i in range(nv):
            for j in range(nv):
                if i == j:
                    eig_v[i, j] = 1.0
                else:
                    den = dft_eig_v[i] - dft_eig_v[j]
                    if den < 1e-12:
                        eig_v[i, j] = 1.0
                    else:
                        num = gw_eig_v[i] - gw_eig_v[j]
                        eig_v[i, j] = num/den 
        
        return eig_c, eig_v 

class Elph:
    def get_elph(self, vbm, nc, nv):
        with h5py.File('./struct_elph_1.h5', 'r') as r:
            elph = np.vectorize(complex)(r['/elph_cart_real'][0, :, 0, : , :], r['/elph_cart_imag'][0, :, 0, : , :]) # g[s\alpha, j, i]
        
        # # Old code. All valence bands are done in this calculation. 
        # elph_c = elph[:, vbm:vbm+nc:1, vbm:vbm+nc:1]
        # elph_v = elph[:, vbm-1:vbm-1-nv:-1, vbm-1:vbm-1-nv:-1]
        
        # New code. elph only includes the bands for absorption calculation. 
        elph_c = elph[:, nv:, nv:]
        elph_v = elph[:, nv-1::-1, nv-1::-1]
        
        return elph_c, elph_v 
    
class Xct:
    def get_xcteig(self):
        with h5py.File('./eigenvectors.h5', 'r') as r:
            xct_eig = r['/exciton_data/eigenvalues'][:][0]
            
        return xct_eig
    
    def get_xctevec(self):
        with h5py.File('./eigenvectors.h5', 'r') as r:
            xct_evec = np.vectorize(complex)(r['/exciton_data/eigenvectors'][0, 0, 0, :, :, 0, 0], r['/exciton_data/eigenvectors'][0, 0, 0, :, :, 0, 1])
            
        return xct_evec  # A[c, v]
    
    def get_vbm_nc_nv(self):
        with h5py.File('./eigenvectors.h5', 'r') as r:
            vbm = r['/mf_header/kpoints/ifmax'][0, 0]
            nc = r['/exciton_header/params/nc'][()]
            nv = r['/exciton_header/params/nv'][()]
            
        return vbm, nc, nv # 1 index based. 
    
class Esf:
    def __init__(self):
        '''
        Class object variables declared and defined in self.assemble_components() function. 
        '''
        pass 
    
    def assemble_components(self):
        self.atoms = DftAtoms().get_atoms()
        
        self.dft_force = DftForce().get_dftforce()
        
        self.xct = Xct()
        self.vbm, self.nc, self.nv = self.xct.get_vbm_nc_nv()
        
        self.eig_c, self.eig_v = Eqp().get_eig(self.vbm, self.nc, self.nv)
        
        self.elph_c, self.elph_v = Elph().get_elph(self.vbm, self.nc, self.nv)
        
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
        self.atoms.calc = EsfCalculator(excited_force=self.excited_force)
        self.atoms.get_potential_energy()
        self.atoms.get_forces()
        with open('esf.xsf', 'w') as w: write_xsf(w, [self.atoms])
        
#endregion

#region: Main.
'''
Only stuff for Gamma point supercell calculations. 
'''
if __name__ == '__main__':
    main()
#endregion
