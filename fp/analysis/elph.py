#region: Modules.
import h5py
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
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
#endregion
