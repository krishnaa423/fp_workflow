#region: Modules.
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EqpResult:
    def __init__(self):
        self.dft: np.ndarray = None  
        self.gw: np.ndarray = None 
    
    def get_dfteig(self):
        data = np.loadtxt('./bandstructure_absorption.dat', skiprows=2, dtype='f8')
        
        self.bands = data[:, 1].astype(dtype='i4')
        self.dft=  data[:, 5]
        
        return self.dft
    
    def get_gweig(self):
        data = np.loadtxt('./bandstructure_absorption.dat', skiprows=2, dtype='f8')
        
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
#endregion
