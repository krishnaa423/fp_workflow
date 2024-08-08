#region: Modules.
import h5py 
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctResult:
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
#endregion
