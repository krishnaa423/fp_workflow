#region: Modules.
import numpy as np 
import matplotlib.pyplot as plt 
import xml.etree.ElementTree as ET
import re 
import glob 
from fp.flows.fullgridflow import FullGridFlow
from fp.io.pkl import load_obj
from fp.structure.kpath import KPath
from scipy.interpolate import griddata
from typing import List, Dict, Tuple
import h5py 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctphPlot:
    def __init__(
        self,
        xctph_filename: str,
        phbands_filename: str,
        fullgridflow_filename: str,
        bandpatkpkl_filename: str,
        xctph_mult_factor: float=1.0,
        xct_state: int=0,   # 0 based index. 
        xct_Qpt_idx: int=0, # 0 based index. 
    ):
        '''
        Inputs:
          xct_state: int
            A zero based index for the exciton state. Default value is 0, which indicates the lowest exciton state. 
        '''
        self.xctph_filename: str = xctph_filename
        self.phbands_filename: str = phbands_filename
        self.fullgridflow_filename: str = fullgridflow_filename
        self.bandpathpkl_filename: str = bandpatkpkl_filename
        self.xctph_mult_factor: float = xctph_mult_factor
        self.xct_state: int = xct_state
        self.xct_Qpt_idx: int = xct_Qpt_idx

        # Additional data created. 
        self.num_bands: int = None 
        self.phbands: np.ndarray = None 
        self.fullgridflow: FullGridFlow = None 
        self.kpath: KPath = None 
        
        self.xctph_interpolated: np.ndarray = None 

    def get_phbands_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[1]
        self.kpath = load_obj(self.bandpathpkl_filename)
        self.fullgridflow = load_obj(self.fullgridflow_filename)

    def get_xctph_data(self):
        ryau2eva = 13.6057039763/0.529177
        xctph: np.ndarray = None 
        qpts: np.ndarray = None 
        with h5py.File(self.xctph_filename, 'r') as r:
            xctph = np.abs(r['xctph'][
                self.xct_state,
                self.xct_state,
                self.xct_Qpt_idx,
                :,
                :
            ]*ryau2eva)
            qpts = r['qpts'][:]
        
        kpath_pts = self.kpath.get_kpts()

        num_kpath_pts = kpath_pts.shape[0]
        num_modes = xctph.shape[0]
        self.xctph_interpolated = np.zeros(shape=(num_kpath_pts, num_modes)) 
        for mode in range(num_modes):
            self.xctph_interpolated[:, mode] = griddata(qpts, xctph[mode, :], kpath_pts, method='linear')*self.xctph_mult_factor

    def save_plot(self, save_filename, show=False, ylim=None):
        # Get some data. 
        self.get_phbands_data()
        self.get_xctph_data()
        path_special_points = self.fullgridflow.path_special_points
        path_segment_npoints = self.fullgridflow.path_segment_npoints

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments or total npoints. 
        if path_segment_npoints:
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
        else:
            xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()    
            ax.plot(xaxis, self.phbands, color='blue')
            xaxis = xaxis.reshape(-1, 1)
            num_modes = self.phbands.shape[1]
            xaxis = np.repeat(xaxis, num_modes, axis=1)
            ax.scatter(xaxis, self.phbands, s=self.xctph_interpolated, color='red')
            ax.yaxis.grid(False) 
            ax.set_xticks(
                ticks=special_points,
                labels=special_labels,
            )

        # Set some labels. 
        ax.set_title(f'Phonon bands and xctph coupling for xct={self.xct_state+1} and Qpt={self.xct_Qpt_idx+1}')
        ax.set_ylabel('Freq (cm-1)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.savefig(save_filename)
        if show: plt.show()
#endregion
