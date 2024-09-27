#region: Modules.
from ase import Atoms 
from ase.dft.kpoints import BandPath, get_special_points
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class KPath:
    def __init__(
            self, 
            atoms: Atoms,
            path_special_points: str, 
            path_segment_npoints: int=None,
            path_total_npoints: int=None,
        ):

        assert (path_segment_npoints is not None) ^ (path_total_npoints is not None) , "Only either path_segments_npoints or path_total_npoints should be present"

        self.atoms: Atoms = atoms
        self.path_special_points: list = path_special_points
        self.path_segment_npoints: int = path_segment_npoints
        self.path_total_npoints: int = path_total_npoints

        # generate bandpath.
        if self.path_total_npoints: 
            self.bandpath: BandPath = atoms.cell.bandpath(path=''.join(self.path_special_points), npoints=len(self.path_special_points)*self.path_segment_npoints)
        else:
            self.bandpath: BandPath = atoms.cell.bandpath(path=''.join(self.path_special_points), npoints=self.path_total_npoints)

    def get_kpts(self):
        
        kpts: np.ndarray = None 
        if self.path_segment_npoints is not None:
            special_points_loc = get_special_points(self.atoms.cell)

            num_special_points = len(self.path_special_points)
            kpts = np.zeros(shape=((num_special_points-1)*self.path_segment_npoints+1, 3), dtype='f8')

            # Add points between the special points. 
            for sp_idx in range(num_special_points-1):
                for coord in range(3):
                    start = special_points_loc[self.path_special_points[sp_idx]][coord]
                    stop = special_points_loc[self.path_special_points[sp_idx+1]][coord]
                    step = (stop - start)/self.path_segment_npoints
                    kpts[sp_idx*self.path_segment_npoints:(sp_idx+1)*self.path_segment_npoints, coord] = np.arange(start, stop, step) if step!=0.0 else 0.0

            # Add the final kpoint. 
            kpts[-1, :] = np.array(special_points_loc[self.path_special_points[-1]])
        else: # path_total_npoints. 
            kpts = self.bandpath.kpts

        return kpts
    
    def get_axis(self):
        
        if self.path_total_npoints is not None:
            return self.bandpath.get_linear_kpoint_axis()
        else: # self.path_segment_npoints
            raise NotImplementedError('Need to implement axis for segments')
    
    def find_K_from_k(self, k, M):
        """Gets a k vector in scaled coordinates and returns a K vector and the
        unfolding G in scaled Coordinates."""

        KG = np.dot(M, k)
        G = np.zeros(3, dtype=int)

        for i in range(3):
            if KG[i] > 0.5000001:
                G[i] = int(np.round(KG[i]))
                KG[i] -= np.round(KG[i])
            elif KG[i] < -0.4999999:
                G[i] = int(np.round(KG[i]))
                KG[i] += abs(np.round(KG[i]))

        return KG, G

    def get_sc_path(self, sc_grid: np.ndarray):
        M = np.diag(sc_grid)
        kpts = self.bandpath.kpts

        Kpts = np.zeros_like(kpts)
        Gpts = np.zeros_like(kpts)
        for kpt_idx, kpt in enumerate(kpts):
            Kpt, G = self.find_K_from_k(kpt, M)
            Kpts[kpt_idx, :] = Kpt
            Gpts[kpt_idx, :] = G

        return Kpts, Gpts
#endregion