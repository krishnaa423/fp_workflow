#region: Modules.
import numpy as np 
import matplotlib.pyplot as plt 
import xml.etree.ElementTree as ET
import re 
import glob 
from fp.inputs.input_main import Input
from fp.io.pkl import load_obj
from fp.structure.kpath import KPath
from scipy.interpolate import griddata
from typing import List, Dict, Tuple
import h5py 
import copy
from scipy import spatial
from ase.units import Ry 
#endregion

#region: Variables.
#endregion

#region: Functions.
def unit_range(klist_in):

    tol = 1e-6
    klist = copy.copy(klist_in)

    for ik, k in enumerate(klist):
        for i, kx in enumerate(k):

            while kx < -tol:
                kx = kx + 1.0
            while kx > 1.0 -tol:
                kx = kx - 1.0

            k[i] = kx

        klist[ik,:] = k

    return klist


def find_kpt(ktargets, klist):

    klist = unit_range(klist)
    ktargets = unit_range(ktargets)

    tree = spatial.KDTree(klist)
    ik_addr = list()
    for k in ktargets:
        d, i = tree.query(k)
        if d > 1e-6:
            print('kpt not found:', k)
            i = None

        ik_addr.append(i)

    return ik_addr
#endregion

#region: Classes.
class XctphPlot:
    def __init__(
        self,
        xctph_filename: str = 'xctph.h5',
        phbands_filename: str = 'struct.freq.gp',
        bandpathpkl_filename: str = 'bandpath.pkl',
        input_filename: str = 'fullgridflow.pkl',
        xctph_mult_factor: float=1.0,
        xct_Qpt_idx: int=0, # 0 based index. 
        xct_state: int=0,   # 0 based index. 
    ):
        '''
        Inputs:
          xct_state: int
            A zero based index for the exciton state. Default value is 0, which indicates the lowest exciton state. 
        '''
        self.xctph_filename: str = xctph_filename
        self.phbands_filename: str = phbands_filename
        self.bandpathpkl_filename: str = bandpathpkl_filename
        self.input_filename: str = input_filename
        self.xctph_mult_factor: float = xctph_mult_factor
        self.xct_state: int = xct_state
        self.xct_Qpt_idx: int = xct_Qpt_idx

        # Additional data created. 
        self.num_bands: int = None 
        self.phbands: np.ndarray = None 
        self.kpath: KPath = None 
        self.input: Input = load_obj(self.input_filename) 
        self.input_dict: dict = self.input.input_dict
        self.xctph_interpolated: np.ndarray = None 

    def get_phbands_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[1]
        self.kpath = load_obj(self.bandpathpkl_filename)

    def get_xctph_gridded(self, xctph_values, kpts_flat):
        kgrid_size = np.array(self.input_dict['wfn']['kdim'], dtype='i4')
        xctph_gridded = np.zeros(shape=kgrid_size)
        for x_idx in range(kgrid_size[0]):
            for y_idx in range(kgrid_size[1]):
                for z_idx in range(kgrid_size[2]):
                    kpt = np.array([
                        x_idx/kgrid_size[0],
                        y_idx/kgrid_size[1],
                        z_idx/kgrid_size[2],
                    ]).reshape(1, 3)
                    kpt_idx = find_kpt(kpt, kpts_flat)[0]
                    if kpt_idx is None or kpt_idx<0:
                        raise Exception(f'kpt_idx is not valid: {kpt_idx}')
                    xctph_gridded[x_idx, y_idx, z_idx] = xctph_values[kpt_idx]

        xctph_gridded = np.pad(xctph_gridded, pad_width=1, mode='wrap')[1:, 1:, 1:]
        kpts_gridded = np.zeros(shape=((kgrid_size[0]+1)*(kgrid_size[1]+1)*(kgrid_size[2]+1), 3))
        x, y, z = np.meshgrid(
            np.linspace(0, 1, kgrid_size[0]+1),
            np.linspace(0, 1, kgrid_size[1]+1),
            np.linspace(0, 1, kgrid_size[2]+1),
        )
        kpts_gridded[:, 0], kpts_gridded[:, 1], kpts_gridded[:, 2] = x.flatten(), y.flatten(), z.flatten()

        return kpts_gridded.reshape(-1, 3), xctph_gridded.reshape(-1, 1)

    def get_xctph_data(self):
        xctph: np.ndarray = None 
        qpts: np.ndarray = None 
        with h5py.File(self.xctph_filename, 'r') as r:
            xctph = np.abs(r['xctph_eh'][
                self.xct_state,
                self.xct_state,
                self.xct_Qpt_idx,
                :,
                :
            ]*Ry)
            qpts = r['qpts'][:]
        
        kpath_pts = self.kpath.get_kpts()

        num_kpath_pts = kpath_pts.shape[0]
        num_modes = xctph.shape[0]
        self.xctph_interpolated = np.zeros(shape=(num_kpath_pts, num_modes)) 
        for mode in range(num_modes):
            kpts_gridded, xctph_gridded = self.get_xctph_gridded(xctph[mode, :], qpts)
            # self.xctph_interpolated[:, mode] = griddata(qpts, xctph[mode, :], kpath_pts, method='linear')*self.xctph_mult_factor
            self.xctph_interpolated[:, mode] = griddata(kpts_gridded, xctph_gridded, kpath_pts, method='linear').reshape(-1)*self.xctph_mult_factor

    def save_plot(self, save_filename='xctph.png', show=False, ylim=None):
        # Get some data. 
        self.get_phbands_data()
        self.get_xctph_data()
        path_special_points = self.kpath.path_special_points
        path_segment_npoints = self.kpath.path_segment_npoints

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments or total npoints. 
        ax.plot(self.phbands, color='blue')
        xaxis = np.arange(self.phbands.shape[0]).reshape(-1, 1)
        num_modes = self.phbands.shape[1]
        xaxis = np.repeat(xaxis, num_modes, axis=1)
        ax.scatter(xaxis, self.phbands, s=self.xctph_interpolated, color='red')
        ax.yaxis.grid(False)  
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )

        # Set some labels. 
        ax.set_title(f'Phonon bands and xctph coupling for xct={self.xct_state} and Qpt={self.xct_Qpt_idx}')
        ax.set_ylabel('Freq (cm-1)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.savefig(save_filename)
        if show: plt.show()
#endregion
