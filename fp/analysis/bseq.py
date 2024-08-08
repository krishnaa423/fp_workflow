#region: Modules.
import os 
from fp.flows.fullgridflow import FullGridFlow
from fp.io import load_obj
from fp.inputs.input_main import Input
from fp.structure import Kgrid
import numpy as np 
import h5py 
from typing import Tuple
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class BseqResult:
    '''
    Will process and get Omega_{Q, s}
    Will process and get A_{Q, s, v, c, k}
    '''
    def __init__(
        self,
        bseq_foldername,
        inputpkl_filename,
    ):
        self.bseq_foldername = bseq_foldername
        self.inputpkl_filename = inputpkl_filename

        self.input: Input = load_obj(self.inputpkl_filename)
        self.Qpts: np.ndarray = None 
        self.xct_eigs: np.ndarray = None 
        self.xct_evecs: np.ndarray = None 

    def get_Qpts(self):
        kdim = self.input.bseq.Qdim
        self.Qpts = Kgrid(
            atoms=self.input.atoms,
            kdim=kdim,
            qshift=(0.0, 0.0, 0.0),
            is_reduced=False,
        ).get_kpts()
        

    def get_foldername_from_Qpt(self, Qpt: np.ndarray) -> str:
        Qpt0 = f'{Qpt[0]:15.10f}'.strip()
        Qpt1 = f'{Qpt[1]:15.10f}'.strip()
        Qpt2 = f'{Qpt[2]:15.10f}'.strip()
        dir_name = f'Q_{Qpt0}_{Qpt1}_{Qpt2}'
        
        return dir_name
    
    def get_Qpt_from_foldername(self, foldername: str) -> np.ndarray:
        Qpt_list = [float(item) for item in foldername.split('_')[1:]]
        
        assert len(Qpt_list)==3, 'Length of the Qpts should be 3 while reading results in bseq folder.'

        Qpt = np.ndarray(Qpt_list).astype('i4')
        return Qpt 

    def get_xct_eigs_and_evecs(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns xct_eigs and xct_evecs. 
        '''
        self.get_Qpts()
        for Qpt_idx, Qpt in enumerate(self.Qpts):
            Qpt_foldername = f'{self.bseq_foldername}/{self.get_foldername_from_Qpt(Qpt)}'
            with h5py.File(f'{Qpt_foldername}/eigenvectors.h5', 'r') as r:
                ds_xct_eigs = r['/exciton_data/eigenvalues']
                ds_xct_evecs = r['/exciton_data/eigenvectors']
                
                # if its the first one, just set the sizes accordingly. 
                if Qpt_idx==0:
                    self.xct_eigs = np.zeros(
                        shape=(
                            self.Qpts.shape[0], # Q. 
                            ds_xct_eigs.shape[0] # x. 
                        ),
                        dtype='f8',
                    ) # Q, s. 
                    self.xct_evecs = np.zeros(
                        shape=(
                            self.Qpts.shape[0],  # Q. 
                            ds_xct_evecs.shape[1], # x
                            ds_xct_evecs.shape[2], # k
                            ds_xct_evecs.shape[3], # c
                            ds_xct_evecs.shape[4], # v
                            ds_xct_evecs.shape[5], # s
                        ),
                        dtype='c16',
                    ) # Q, x, k, c, v, s. 
            
                # Get the datasets. 
                self.xct_eigs[Qpt_idx, :] = ds_xct_eigs[:]
                self.xct_evecs[Qpt_idx, ...] = np.vectorize(complex)(ds_xct_evecs[0, :, :, :, :, :, 0], ds_xct_evecs[0, :, :, :, :, :, 1])

        return self.xct_eigs, self.xct_evecs

    def bseq_results(self):
        self.get_xct_eigs_and_evecs()
#endregion
