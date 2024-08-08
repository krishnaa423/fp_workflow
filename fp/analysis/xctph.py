#region: Modules.
from fp.analysis import *
import numpy as np 
from scipy.spatial import KDTree
import h5py 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPh:
    def __init__(
        self,
        elph_files_prefix,
        bseq_foldername, 
        fullgridflowpkl_filename,
        inputpkl_filename,
    ):
        self.elph_files_prefix = elph_files_prefix
        self.bseq_foldername = bseq_foldername
        self.fullgridflowpkl_filename = fullgridflowpkl_filename
        self.inputpkl_filename = inputpkl_filename

        # Will be populated by class functions. 
        self.elph_c: np.ndarray = None 
        self.elph_v: np.ndarray = None 
        self.ph_eigs: np.ndarray = None 
        self.xct_eigs: np.ndarray = None 
        self.xct_evecs: np.ndarray = None 
        self.kpts: np.ndarray = None 
        self.qpts: np.ndarray = None 
        self.Qpts: np.ndarray = None 
        self.num_kpts: int = None 
        self.kk_map: np.ndarray = None 
        self.xct_evecs_Qq: np.ndarray = None 
        self.xct_evecs_kq: np.ndarray = None 
        self.elph_c_kQ: np.ndarray = None 
        self.xctph: np.ndarray = None 

    def get_add_idxs(self, a, b, kdtree: KDTree) -> np.ndarray:
        frac, _ = np.modf(a + b)
        _, idxs = kdtree.query(frac)
        return idxs
        
    def get_kpts_plus_kpts_map(self):
        self.kk_map = np.zeros(shape=(self.num_kpts, self.num_kpts), dtype='i4')
        kdtree = KDTree(self.kpts)
        k1 = np.repeat(self.kpts, self.num_kpts, axis=0)
        k2 = np.tile(self.kpts, reps=(self.num_kpts, 1))
        idxs = self.get_add_idxs(k1, k2, kdtree=kdtree)
        self.kk_map = idxs.reshape(self.num_kpts, self.num_kpts)

        return self.kk_map

    def get_xct_evecs_Qq(self):
        self.xct_evecs_Qq = np.zeros(shape=(*self.xct_evecs.shape, self.num_kpts), dtype='c16')
        for Qidx in range(self.num_kpts):
            for qidx in range(self.num_kpts):
                self.xct_evecs_Qq[Qidx, :, :, :, :, :, qidx] = self.xct_evecs[self.kk_map[Qidx, qidx], ...]

    def get_xct_evecs_kq(self):
        self.xct_evecs_kq = np.zeros(shape=(*self.xct_evecs.shape, self.num_kpts), dtype='c16')
        for kidx in range(self.num_kpts):
            for qidx in range(self.num_kpts):
                self.xct_evecs_kq[:, :, kidx, :, :, :, qidx] = self.xct_evecs[:, :, self.kk_map[kidx, qidx], :, :, :]

    def get_elph_c_kQ(self):
        self.elph_c_kQ = np.zeros(shape=(*self.elph_c.shape, self.num_kpts), dtype='c16')
        for kidx in range(self.num_kpts):
            for Qidx in range(self.num_kpts):
                self.elph_c_kQ[:, :, kidx, :, :, Qidx] = self.elph_c[:, :, self.kk_map[kidx, Qidx], :, :]

    def assemble_components(self):
        # xct_eig, xct_evec. 
        bseq_result = BseqResult(
            bseq_foldername=self.bseq_foldername,
            inputpkl_filename=self.inputpkl_filename,
        )
        self.xct_eigs, self.xct_evecs = bseq_result.get_xct_eigs_and_evecs()
        

        # elph_c, elph_v.
        elphq_result = ElphQResult(
            elph_files_prefix=self.elph_files_prefix,
            fullgridflowpkl_filename=self.fullgridflowpkl_filename,
        )
        self.elph_c, self.elph_v = elphq_result.get_elph(ev_units=False)
        self.ph_eigs = elphq_result.ph_eigs

        # kpts, qpts, Qpts. 
        kpts_result = KptsResult(
            self.inputpkl_filename
        )
        self.kpts = kpts_result.get_wfn_kpts()
        self.num_kpts = self.kpts.shape[0]

        # All the idxs. 
        self.get_kpts_plus_kpts_map()

        # Create additional add tensors. 
        self.get_xct_evecs_kq()
        self.get_xct_evecs_Qq()
        self.get_elph_c_kQ()

    def get_xctph(self, zero_el=False):
        self.assemble_components()

        num_Qpts = self.xct_eigs.shape[0]
        num_xct_eigs = self.xct_eigs.shape[1]
        num_phmodes = self.elph_c.shape[1]
        # num_qpts = self.Qpts      # Assume same. 

        self.xctph = np.zeros(shape=(num_xct_eigs, num_xct_eigs, num_phmodes, num_Qpts, num_Qpts), dtype='c16')

        # Parts from eq (9) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        if not zero_el:
            self.xctph += np.einsum(
                'axkcvsq,qukCca,aXkCvs->xXuaq',
                np.conjugate(self.xct_evecs_Qq),
                self.elph_c_kQ,
                self.xct_evecs,
            )

        self.xctph += np.einsum(
            'axkcvsq,qukVv,aXkcVsq->xXuaq',
            np.conjugate(self.xct_evecs_Qq),
            self.elph_v,
            self.xct_evecs_kq,
        )

        return self.xctph

    def write(self):

        with h5py.File('xctph.h5', 'w') as w:
            ds_xctph = w.create_dataset('xctph', shape=self.xctph.shape, dtype=self.xctph.dtype)
            ds_xctph[:] = self.xctph

#endregion
