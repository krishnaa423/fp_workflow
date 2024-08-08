#region: Modules.
import h5py
import numpy as np 
import glob
from fp.io import load_obj
from fp.flows.fullgridflow import FullGridFlow
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class ElphQResult:
    def __init__(
        self,
        elph_files_prefix,   
        fullgridflowpkl_filename,
    ):
        self.elph_files_prefix: str = elph_files_prefix
        self.fullgridflowpkl_filename: str = fullgridflowpkl_filename

        # Init processing. 
        self.fullgridflow: FullGridFlow = load_obj(self.fullgridflowpkl_filename)

        # During run. 
        self.elph_files: list = None 
        self.elph_c: np.ndarray = None 
        self.elph_v: np.ndarray = None 
        self.ph_eigs: np.ndarray = None 

    def get_elph_files(self):
        # Assumption is after sorting it would have suffixes 1, 2, 3... and so on. 
        self.elph_files = glob.glob(f'{self.elph_files_prefix}*')
        self.elph_files.sort()

        num_kpts = int(np.prod(self.fullgridflow.wfn_qe_kdim))

        assert num_kpts%len(self.elph_files)==0, f'Number of kpoints {num_kpts} does not evenly divide among number of elph files {len(self.elph_files)}.'
        
        return self.elph_files

    def get_elph(self, ev_units=False):
        '''
        Get elph_c and elph_v. 
        '''
        elph_files = self.get_elph_files()
        num_elph_files = len(elph_files)
        num_kpts = int(np.prod(self.fullgridflow.wfn_qe_kdim))
        num_kpts_per_elph_file = int(num_kpts//num_elph_files)
        num_cond = self.fullgridflow.abs_cond_bands
        num_val = self.fullgridflow.abs_val_bands

        #region: To get sizes initialized. 
        ds_elph_shape = None 
        with h5py.File(elph_files[0], 'r') as r:
            ds_elph_sample = r['/elph_cart_real']
            ds_elph_shape = ds_elph_sample.shape

        self.elph_c = np.zeros(
            shape=(
                num_kpts, # q
                ds_elph_shape[1], # nu
                num_kpts, # k
                num_cond, # c'
                num_cond, # c
            ),
            dtype='c16',
        )

        self.elph_v = np.zeros(
            shape=(
                num_kpts, # q
                ds_elph_shape[1], # nu
                num_kpts, # k
                num_val, # v'
                num_val, # v
            ),
            dtype='c16',
        )
        #endregion.

        # Populate the arrays now. 
        for elph_file_idx in range(num_elph_files):
            with h5py.File(elph_files[elph_file_idx], 'r') as r:
                elph = np.vectorize(complex)(r['/elph_nu_real'][:], r['/elph_nu_imag'][:])       # g[q, nu, k, j, i]
                if ev_units: elph *= 13.605693122994 # Ry -> eV. 
                self.elph_c[:] = elph[:, :, elph_file_idx*num_kpts_per_elph_file:(elph_file_idx+1)*num_kpts_per_elph_file, num_val:, num_val:]
                self.elph_v[:] = elph[:, :, elph_file_idx*num_kpts_per_elph_file:(elph_file_idx+1)*num_kpts_per_elph_file, num_val-1::-1, num_val-1::-1]

        # Get ph_eigs. 
        with h5py.File(self.elph_files[0], 'r') as r:
            self.ph_eigs = r['/ph_eigs'][:]

        return self.elph_c, self.elph_v

class ElphResult:
    def get_elph(self, vbm, nc, nv):
        with h5py.File('./struct_elph_1.h5', 'r') as r:
            elph = np.vectorize(complex)(r['/elph_cart_real'][0, :, 0, : , :], r['/elph_cart_imag'][0, :, 0, : , :])*27.2114/0.529177      # Ha/bohr -> eV/A.  # g[s\alpha, j, i]
        
        # # Old code. All valence bands are done in this calculation. 
        # elph_c = elph[:, vbm:vbm+nc:1, vbm:vbm+nc:1]
        # elph_v = elph[:, vbm-1:vbm-1-nv:-1, vbm-1:vbm-1-nv:-1]
        
        # New code. elph only includes the bands for absorption calculation. 
        elph_c = elph[:, nv:, nv:]
        elph_v = elph[:, nv-1::-1, nv-1::-1]
        
        return elph_c, elph_v 
#endregion
