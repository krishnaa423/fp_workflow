#region: Modules.
import numpy as np 
import h5py 
from scipy.spatial import KDTree
from xctph.kpoints import get_all_kq_maps
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPolCalc:
    def __init__(
        self,
        eph_filename,
        xctph_filename,
        xctph_partial_filename,
        fullgridflow_filename, 
        max_error=1.0e-4,
        max_steps=10,
    ):
        self.eph_filename = eph_filename
        self.xctph_filename = xctph_filename
        self.xctph_partial_filename = xctph_partial_filename
        self.fullgridflow_filename = fullgridflow_filename
        self.max_error = max_error
        self.max_steps = max_steps 

        # Later. 
        self.xctph_zeroel: np.ndarray = None 
        self.xctph: np.ndarray = None 
        self.xct_eigs: np.ndarray = None 
        self.ph_eigs: np.ndarray = None 
        self.kpts: np.ndarray = None 
        self.num_kpts: int = None 
        self.num_xct_eigs: int = None 
        self.num_phmodes: int = None 

        self.kk_plus_map: np.ndarray = None 
        self.kk_minus_map: np.ndarray = None 
        self.xctph_QQ: np.ndarray = None 
        self.xctpol_el_QQ: np.ndarray = None 
        self.xctpol_ph_QQ: np.ndarray = None 

        self.xctpol_el_eig_previous: float = None  # A in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        self.xctpol_el_eig: float = None  # A in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        self.xctpol_el: np.ndarray = None  # A in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        self.xctpol_ph: np.ndarray = None  # B in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        self.xctpol_energy: float = None 
        self.eigs: np.ndarray = None 
        self.evecs: np.ndarray = None 

    def get_add_idxs(self, a, b, kdtree: KDTree) -> np.ndarray:
        frac, _ = np.modf(a + b)
        _, idxs = kdtree.query(frac)
        return idxs
    
    def get_sub_idxs(self, a, b, kdtree: KDTree) -> np.ndarray:
        sub_ab = a - b
        sub_ab[sub_ab<0.0] += 1.0
        frac, _ = np.modf(sub_ab)
        _, idxs = kdtree.query(frac)
        return idxs
        
    def get_kk_plus_map(self):
        # self.kk_plus_map = np.zeros(shape=(self.num_kpts, self.num_kpts), dtype='i4')
        # kdtree = KDTree(self.kpts)
        # k1 = np.repeat(self.kpts, self.num_kpts, axis=0)
        # k2 = np.tile(self.kpts, reps=(self.num_kpts, 1))
        # idxs = self.get_add_idxs(k1, k2, kdtree=kdtree)
        # self.kk_plus_map = idxs.reshape(self.num_kpts, self.num_kpts)

        # Using algo from xctph package. 
        self.kk_plus_map = get_all_kq_maps(self.kpts, self.kpts, plus_or_minus=1.0)

        return self.kk_plus_map

    def get_kk_minus_map(self):
        # self.kk_minus_map = np.zeros(shape=(self.num_kpts, self.num_kpts), dtype='i4')
        # kdtree = KDTree(self.kpts)
        # k1 = np.repeat(self.kpts, self.num_kpts, axis=0)
        # k2 = np.tile(self.kpts, reps=(self.num_kpts, 1))
        # idxs = self.get_sub_idxs(k1, k2, kdtree=kdtree)
        # self.kk_minus_map = idxs.reshape(self.num_kpts, self.num_kpts)

        # Using algo from xctph package. 
        self.kk_minus_map = get_all_kq_maps(self.kpts, self.kpts, plus_or_minus=-1.0)

        return self.kk_minus_map
    
    def assemble_components(self):


        with h5py.File(self.xctph_partial_filename, 'r') as r:
            self.xctph_zeroel = r['/xctph'][:]      # i, j, Q, nu, q. 
            self.xctph_zeroel = np.einsum('ijauq->ijuaq', self.xctph_zeroel) # i, j, nu, Q, q. 

        with h5py.File(self.xctph_filename, 'r') as r:
            self.xctph = r['/xctph'][:]      # i, j, Q, nu, q. 
            self.xctph = np.einsum('ijauq->ijuaq', self.xctph) # i, j, nu, Q, q. 
            
            self.xct_eigs = r['/energies'][:] # S, Q.  
            self.xct_eigs = np.einsum('SQ->QS', self.xct_eigs)  # Q, S. 

            self.kpts = r['/Qpts'][:] # nk, 3.
            self.num_kpts = self.kpts.shape[0]
            self.num_xct_eigs = self.xct_eigs.shape[1]
        
        with h5py.File(self.eph_filename, 'r') as r:
            self.ph_eigs = r['/gkq_data/frequencies'][:] # nu, q 
            self.ph_eigs = np.einsum('uq->qu', self.ph_eigs) # q, nu. 

            self.num_phmodes = self.ph_eigs.shape[1]

        # get maps. 
        self.get_kk_plus_map()
        self.get_kk_minus_map()

        # Threshold the ph_eigs array. Just doing something random. 
        self.ph_eigs[np.abs(self.ph_eigs) < 1e-12] = 1e-6

    def get_xctpol_el_QQ(self):
        self.xctpol_el_QQ = np.zeros(shape=(self.num_kpts, self.num_xct_eigs ,self.num_kpts), dtype='c16')
        for Qidx in range(self.num_kpts):
            for Qpidx in range(self.num_kpts):
                self.xctpol_el_QQ[Qidx, :, Qpidx] = self.xctpol_el[self.kk_plus_map[Qidx, Qpidx], :]

    def get_xctpol_ph_QQ(self):
        self.xctpol_ph_QQ = np.zeros(shape=(self.num_kpts, self.num_phmodes ,self.num_kpts), dtype='c16')
        for Qidx in range(self.num_kpts):
            for Qpidx in range(self.num_kpts):
                self.xctpol_ph_QQ[Qidx, :, Qpidx] = self.xctpol_ph[self.kk_minus_map[Qidx, Qpidx], :]

    def get_xctph_QQ(self):
        self.xctph_QQ = np.zeros(shape=(self.num_xct_eigs, self.num_xct_eigs, self.num_phmodes, self.num_kpts, self.num_kpts, self.num_kpts), dtype='c16')
        for Qidx in range(self.num_kpts):
            for Qpidx in range(self.num_kpts):
                self.xctph_QQ[:, :, :, :, Qidx, Qpidx] = self.xctph[:, :, :, :, self.kk_minus_map[Qidx, Qpidx]]

    def solve_first_equation(self):
        '''
        Is an eigensystem. Refer to eq (7) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        '''
        matrix: np.ndarray = np.zeros(shape=(self.num_kpts, self.num_xct_eigs, self.num_kpts, self.num_xct_eigs), dtype='c16')
        
        # Assemble the diagonal part. 
        delta_xct_eigs = np.eye(self.num_xct_eigs)
        delta_kpts = np.eye(self.num_kpts)
        matrix += np.einsum(
            'ax,xX,aA->axAX',
            self.xct_eigs,
            delta_xct_eigs,
            delta_kpts,
        )

        # Assemble the other part.
        self.get_xctpol_ph_QQ()
        self.get_xctph_QQ() 
        matrix -= 2*np.einsum(
            'auA,xXuAaA->axAX',
            self.xctpol_ph_QQ,
            self.xctph_QQ,
        )

        matrix = matrix.reshape(self.num_kpts*self.num_xct_eigs, self.num_kpts*self.num_xct_eigs)

        self.eigs, self.evecs = np.linalg.eig(matrix)
        min_index = np.argmin(self.eigs)
        self.xctpol_el_eig_previous = self.xctpol_el_eig
        self.xctpol_el_eig = self.eigs[min_index].real
        self.xctpol_el = self.evecs[:, min_index].reshape(self.num_kpts, self.num_xct_eigs)

    def solve_second_equation(self, use_zeroel=False):
        '''
        Just an evaluation. Refer to eq (8) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        '''
        # Get xctpol_ph. 
        self.get_xctpol_el_QQ()
        xctpol_ph_uptofactor = np.einsum(
            'AX,axA,xXuAa->au',
            np.conjugate(self.xctpol_el),
            self.xctpol_el_QQ,
            np.conjugate(self.xctph_zeroel if use_zeroel else self.xctph),
        )
        # self.xctpol_ph = xctpol_ph_uptofactor/self.ph_eigs/self.num_kpts
        self.xctpol_ph = xctpol_ph_uptofactor/self.ph_eigs

    def initial_guess(self):
        '''
        Obtained by setting elph_c to zero. 
        '''
        initial_constant = np.sqrt(1/self.num_xct_eigs) # Based on page 12 para below fig 5. 
        self.xctpol_el = np.ones(shape=(self.num_kpts, self.num_xct_eigs), dtype='c16')*initial_constant

        self.solve_second_equation(use_zeroel=True)

        self.solve_first_equation()

    def get_error(self):
        return abs(self.xctpol_el_eig - self.xctpol_el_eig_previous)

    def run(self):

        self.solve_second_equation()
        self.solve_first_equation()

        error = self.get_error()
        step_idx = 0

        # Loop till convergence. 
        while error>self.max_error and step_idx<self.max_steps:
            print(f'Step: {step_idx}, Solving second equation', flush=True)
            self.solve_second_equation()
            print(f'Step: {step_idx}, Done Solving second equation', flush=True)
            
            print(f'Step: {step_idx}, Solving first equation', flush=True)
            self.solve_first_equation()
            print(f'Step: {step_idx}, Done Solving first equation', flush=True)
            
            print(f'Step: {step_idx}, Getting error', flush=True)
            error = self.get_error()
            print(f'Step: {step_idx}, Done getting error. Error: {error}', flush=True)

            print('\n\n\n\n', flush=True)

            step_idx += 1

        if step_idx>=self.max_steps:
            raise Exception(f'Max steps: {step_idx} reached before convergence. Error: {error}')

    #region: Main user interface. 
    def get_xctpol(self):
        print('Assembling components', flush=True)
        self.assemble_components()
        print('Done assembling components', flush=True)
        
        print('Getting initial guess', flush=True)
        self.initial_guess()
        print('Done getting initial guess', flush=True)
        
        print('Running', flush=True)
        self.run()
        print('Done with xctpol', flush=True)

    def get_xctpol_energy(self, evUnits=True):
        '''
        From eq (10) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
        '''
        print('Calculating xctpol energy', flush=True)
        xctpol_el_sq = np.abs(self.xctpol_el)**2
        xctpol_ph_sq = np.abs(self.xctpol_ph)**2

        xctpol_energy = 0
        xctpol_energy += np.einsum(
            'xa,xa->',
            xctpol_el_sq,
            self.xct_eigs,
        )/self.num_kpts

        xctpol_energy -= np.einsum(
            'qu,qu->',
            xctpol_ph_sq,
            self.ph_eigs   
        )/self.num_kpts

        if evUnits: xctpol_energy *= 13.6057039763

        self.xctpol_energy = xctpol_energy

        print('Done calculating xctpol energy', flush=True)

        return self.xctpol_energy
    
    def write(self):

        ry2ev = 13.6057039763

        print('Writing xctpol.h5', flush=True)
        with h5py.File('xctpol.h5', 'w') as w:
            w.create_dataset('xctpol_el', data=self.xctpol_el)
            w.create_dataset('xctpol_ph', data=self.xctpol_ph)
            w.create_dataset('xctpol_energy', data=self.xctpol_energy)
            w.create_dataset('xctpol_eigs', data=self.eigs*ry2ev)
            w.create_dataset('xctpol_eigs_ordered', data=np.sort(self.eigs.real)*ry2ev)
            w.create_dataset('xctpol_evecs', data=self.evecs)
        print('Done writing xctpol.h5', flush=True)
    #endregion

# region: Old version. 
# class XctPol:
#     def __init__(
#         self,
#         elph_files_prefix,
#         bseq_foldername, 
#         fullgridflowpkl_filename,
#         inputpkl_filename,
#         max_error=1.0e-4,
#         max_steps=10,
#     ):
#         self.elph_files_prefix = elph_files_prefix
#         self.bseq_foldername = bseq_foldername
#         self.fullgridflowpkl_filename = fullgridflowpkl_filename
#         self.inputpkl_filename = inputpkl_filename
#         self.max_error = max_error
#         self.max_steps = max_steps 

#         # Later. 
#         self.xctph_container: XctPh = None 
#         self.xctph_zeroel: np.ndarray = None 
#         self.xctph: np.ndarray = None 
#         self.xct_eigs: np.ndarray = None 
#         self.ph_eigs: np.ndarray = None 
#         self.kpts: np.ndarray = None 
#         self.num_kpts: int = None 
#         self.num_xct_eigs: int = None 
#         self.num_phmodes: int = None 

#         self.kk_plus_map: np.ndarray = None 
#         self.kk_minus_map: np.ndarray = None 
#         self.xctph_QQ: np.ndarray = None 
#         self.xctpol_el_QQ: np.ndarray = None 
#         self.xctpol_ph_QQ: np.ndarray = None 

#         self.xctpol_el_eig_previous: float = None  # A in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         self.xctpol_el_eig: float = None  # A in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         self.xctpol_el: np.ndarray = None  # A in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         self.xctpol_ph: np.ndarray = None  # B in the paper. PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         self.xctpol_energy: float = None 

#     def get_add_idxs(self, a, b, kdtree: KDTree) -> np.ndarray:
#         frac, _ = np.modf(a + b)
#         _, idxs = kdtree.query(frac)
#         return idxs
    
#     def get_sub_idxs(self, a, b, kdtree: KDTree) -> np.ndarray:
#         sub_ab = a - b
#         sub_ab[sub_ab<0.0] += 1.0
#         frac, _ = np.modf(sub_ab)
#         _, idxs = kdtree.query(frac)
#         return idxs
        
#     def get_kk_plus_map(self):
#         self.kk_plus_map = np.zeros(shape=(self.num_kpts, self.num_kpts), dtype='i4')
#         kdtree = KDTree(self.kpts)
#         k1 = np.repeat(self.kpts, self.num_kpts, axis=0)
#         k2 = np.tile(self.kpts, reps=(self.num_kpts, 1))
#         idxs = self.get_add_idxs(k1, k2, kdtree=kdtree)
#         self.kk_plus_map = idxs.reshape(self.num_kpts, self.num_kpts)

#         return self.kk_plus_map

#     def get_kk_minus_map(self):
#         self.kk_minus_map = np.zeros(shape=(self.num_kpts, self.num_kpts), dtype='i4')
#         kdtree = KDTree(self.kpts)
#         k1 = np.repeat(self.kpts, self.num_kpts, axis=0)
#         k2 = np.tile(self.kpts, reps=(self.num_kpts, 1))
#         idxs = self.get_sub_idxs(k1, k2, kdtree=kdtree)
#         self.kk_minus_map = idxs.reshape(self.num_kpts, self.num_kpts)

#         return self.kk_minus_map

#     def assemble_components(self):
#         self.xctph_container: XctPh = XctPh(
#             elph_files_prefix=self.elph_files_prefix,
#             bseq_foldername=self.bseq_foldername,
#             fullgridflowpkl_filename=self.fullgridflowpkl_filename,
#             inputpkl_filename=self.inputpkl_filename,
#         )
#         self.xctph_zeroel = self.xctph_container.get_xctph(zero_el=True)
#         self.xctph = self.xctph_container.get_xctph(zero_el=False)
#         self.xct_eigs = self.xctph_container.xct_eigs
#         self.ph_eigs = self.xctph_container.ph_eigs
#         self.kpts = self.xctph_container.kpts
#         self.num_kpts = self.xctph_container.num_kpts
#         self.num_xct_eigs = self.xctph_container.xct_eigs.shape[1]
#         self.num_phmodes = self.xctph_container.elph_c.shape[1]

#         # get maps. 
#         self.get_kk_plus_map()
#         self.get_kk_minus_map()

#         # Threshold the ph_eigs array. Just doing something random. 
#         self.ph_eigs[np.abs(self.ph_eigs) < 1e-12] = 1e-6

#     def get_xctpol_el_QQ(self):
#         self.xctpol_el_QQ = np.zeros(shape=(self.num_kpts, self.num_xct_eigs ,self.num_kpts), dtype='c16')
#         for Qidx in range(self.num_kpts):
#             for Qpidx in range(self.num_kpts):
#                 self.xctpol_el_QQ[Qpidx, :, Qidx] = self.xctpol_el[self.kk_plus_map[Qpidx, Qidx], :]

#     def get_xctpol_ph_QQ(self):
#         self.xctpol_ph_QQ = np.zeros(shape=(self.num_kpts, self.num_phmodes ,self.num_kpts), dtype='c16')
#         for Qidx in range(self.num_kpts):
#             for Qpidx in range(self.num_kpts):
#                 self.xctpol_ph_QQ[Qidx, :, Qpidx] = self.xctpol_ph[self.kk_minus_map[Qidx, Qpidx], :]

#     def get_xctph_QQ(self):
#         self.xctph_QQ = np.zeros(shape=(self.num_xct_eigs, self.num_xct_eigs, self.num_phmodes, self.num_kpts, self.num_kpts, self.num_kpts), dtype='c16')
#         for Qidx in range(self.num_kpts):
#             for Qpidx in range(self.num_kpts):
#                 self.xctph_QQ[:, :, :, :, Qidx, Qpidx] = self.xctph[:, :, :, :, self.kk_minus_map[Qidx, Qpidx]]

#     def solve_first_equation(self):
#         '''
#         Is an eigensystem. Refer to eq (7) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         '''
#         matrix: np.ndarray = np.zeros(shape=(self.num_kpts, self.num_xct_eigs, self.num_kpts, self.num_xct_eigs), dtype='c16')
        
#         # Assemble the diagonal part. 
#         delta_xct_eigs = np.eye(self.num_xct_eigs)
#         delta_kpts = np.eye(self.num_kpts)
#         matrix += np.einsum(
#             'ax,xX,aA->axAX',
#             self.xct_eigs,
#             delta_xct_eigs,
#             delta_kpts,
#         )

#         # Assemble the other part.
#         self.get_xctpol_ph_QQ()
#         self.get_xctph_QQ() 
#         matrix -= 2/self.num_kpts*np.einsum(
#             'auA,xXuAaA->axAX',
#             self.xctpol_ph_QQ,
#             self.xctph_QQ,
#         )

#         matrix = matrix.reshape(self.num_kpts*self.num_xct_eigs, self.num_kpts*self.num_xct_eigs)

#         eigs, evecs = np.linalg.eig(matrix)
#         min_index = np.argmin(eigs)
#         self.xctpol_el_eig_previous = self.xctpol_el_eig
#         self.xctpol_el_eig = eigs[min_index].real
#         self.xctpol_el = evecs[:, min_index].reshape(self.num_kpts, self.num_xct_eigs)

#     def solve_second_equation(self, use_zeroel=False):
#         '''
#         Just an evaluation. Refer to eq (8) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         '''
#         # Get xctpol_ph. 
#         self.get_xctpol_el_QQ()
#         xctpol_ph_uptofactor = np.einsum(
#             'AX,Axa,xXuAa->au',
#             np.conjugate(self.xctpol_el),
#             self.xctpol_el_QQ,
#             np.conjugate(self.xctph_zeroel if use_zeroel else self.xctph),
#         )
#         self.xctpol_ph = xctpol_ph_uptofactor/self.ph_eigs/self.num_kpts

#     def initial_guess(self):
#         '''
#         Obtained by setting elph_c to zero. 
#         '''
#         initial_constant = np.sqrt(1/self.num_xct_eigs) # Based on page 12 para below fig 5. 
#         self.xctpol_el = np.ones(shape=(self.num_kpts, self.num_xct_eigs), dtype='c16')*initial_constant

#         self.solve_second_equation(use_zeroel=True)

#         self.solve_first_equation()

#     def get_error(self):
#         return abs(self.xctpol_el_eig - self.xctpol_el_eig_previous)

#     def run(self):

#         self.solve_second_equation()
#         self.solve_first_equation()

#         error = self.get_error()
#         step_idx = 0

#         # Loop till convergence. 
#         while error>self.max_error and step_idx<self.max_steps:
#             print(f'Step: {step_idx}, Solving second equation', flush=True)
#             self.solve_second_equation()
#             print(f'Step: {step_idx}, Done Solving second equation', flush=True)
            
#             print(f'Step: {step_idx}, Solving first equation', flush=True)
#             self.solve_first_equation()
#             print(f'Step: {step_idx}, Done Solving first equation', flush=True)
            
#             print(f'Step: {step_idx}, Getting error', flush=True)
#             error = self.get_error()
#             print(f'Step: {step_idx}, Done getting error. Error: {error}', flush=True)

#             print('\n\n\n\n', flush=True)

#             step_idx += 1

#         if step_idx>=self.max_steps:
#             raise Exception(f'Max steps: {step_idx} reached before convergence. Error: {error}')

#     #region: Main user interface. 
#     def get_xctpol(self):
#         print('Assembling components', flush=True)
#         self.assemble_components()
#         print('Done assembling components', flush=True)
        
#         print('Getting initial guess', flush=True)
#         self.initial_guess()
#         print('Done getting initial guess', flush=True)
        
#         print('Running', flush=True)
#         self.run()
#         print('Done with xctpol', flush=True)

#     def get_xctpol_energy(self, evUnits=True):
#         '''
#         From eq (10) in PHYSICAL REVIEW LETTERS 132, 036902 (2024)
#         '''
#         print('Calculating xctpol energy', flush=True)
#         xctpol_el_sq = np.abs(self.xctpol_el)**2
#         xctpol_ph_sq = np.abs(self.xctpol_ph)**2

#         xctpol_energy = 0
#         xctpol_energy += np.einsum(
#             'xa,xa->',
#             xctpol_el_sq,
#             self.xct_eigs,
#         )/self.num_kpts

#         xctpol_energy -= np.einsum(
#             'qu,qu->',
#             xctpol_ph_sq,
#             self.ph_eigs   
#         )/self.num_kpts

#         if evUnits: xctpol_energy *= 13.6057039763

#         self.xctpol_energy = xctpol_energy

#         print('Done calculating xctpol energy', flush=True)

#         return self.xctpol_energy
    
#     def write(self):
#         print('Writing xctpol.h5', flush=True)
#         with h5py.File('xctpol.h5', 'w') as w:
#             ds_xctpol_el = w.create_dataset('xctpol_el', shape=self.xctpol_el.shape, dtype=self.xctpol_el.dtype)
#             ds_xctpol_el[:] = self.xctpol_el

#             ds_xctpol_ph = w.create_dataset('xctpol_ph', shape=self.xctpol_ph.shape, dtype=self.xctpol_ph.dtype)
#             ds_xctpol_ph[:] = self.xctpol_ph

#             ds_xctpol_energy = w.create_dataset('xctpol_energy', data=self.xctpol_energy)
#         print('Done writing xctpol.h5', flush=True)
#     #endregion
# endregion: Old version. 
#endregion
