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
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class KpdosPlot:
    def __init__(
        self,
        dftelbands_xml_filename: str,
        fullgridflow_filename: str,
        bandpathpkl_filename,
        kpdos_prefix: str,
        kpdos_mult_factor: float=1.0,
    ):
        self.dftelbands_xml_filename: str = dftelbands_xml_filename
        self.kpdos_prefix: str = kpdos_prefix
        self.bandpathpkl_filename = bandpathpkl_filename
        self.fullgridflow_filename = fullgridflow_filename
        self.kpdos_mult_factor = kpdos_mult_factor

        # Data.
        self.dft_eigs: np.ndarray = None 
        self.fermi_energy: float = None 
        self.fullgridflow: FullGridFlow = None 
        self.kpath: KPath = None 
        self.is_axis_calculated: bool = False
        
        self.kpts: np.ndarray = None 
        self.energy: np.ndarray = None 
        self.num_kpts: int = None 
        self.num_eidxs: int = None  
        self.K: np.ndarray = None 
        self.E: np.ndarray = None 

    def get_annotated_file_list(self) -> List[Dict]:
        files = glob.glob(self.kpdos_prefix)
        
        annotated_file_list: List[Dict] = []
        for file in files:
            match = re.search(r'.*?_atm#\d+\((?P<atom>.*?)\)_wfc#\d+\((?P<orbital>.*?)\)', file)
            annotated_file_list.append({
                'filename': file,
                'atom': match.group('atom'),
                'orbital': match.group('orbital'),
            })

        return annotated_file_list

    def get_dftelbands_data(self):
        ha2eV = 27.2114
        tree = ET.parse(self.dftelbands_xml_filename)
        root = tree.getroot()

        eig_nodes = root.findall('.//ks_energies/eigenvalues')
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*ha2eV
        num_kpts = len(eig_nodes)
        num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
        dft_eigs = np.zeros(shape=(num_kpts, num_bands), dtype='f8')
        for kpt_idx, node in enumerate(eig_nodes):
            dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*ha2eV - fermi_energy

        self.dft_eigs = dft_eigs
        self.fermi_energy = fermi_energy

        self.kpath = load_obj(self.bandpathpkl_filename)
        self.fullgridflow = load_obj(self.fullgridflow_filename)

    def get_kpdos_data(self, filename) -> Tuple[np.ndarray, np.ndarray]:
        data = np.loadtxt(filename, skiprows=1) 

        if not self.is_axis_calculated:
            self.kpts, kpt_counts = np.unique(data[:, 0], return_counts=True)
            self.num_kpts = self.kpts.size
            self.num_eidxs = kpt_counts[0]

            self.kpts = np.arange(self.num_kpts)
            self.energy = data[:, 0:self.num_eidxs] - self.fermi_energy # in eV.
            self.K, self.E = np.meshgrid(self.kpts, self.energy, indexing='ij')
        
        pdos = data[:, -1].reshape(self.num_kpts, -1)*self.kpdos_mult_factor
        pdos = griddata(
            points=np.stack((self.K.reshape(-1, 1), self.E.reshape(-1, 1)), axis=1),
            values=pdos.flatten(),
            xi=np.stack((self.K.reshape(-1, 1), self.dft_eigs.reshape(-1, 1)), axis=1),
            method='linear',
        ).reshape(self.num_kpts, -1)

        return pdos

    def save_plot(self, save_filename, show=False, ylim=None):
        self.get_dftelbands_data()
        annotated_file_list = self.get_annotated_file_list()
        path_special_points = self.fullgridflow.path_special_points
        path_segment_npoints = self.fullgridflow.path_segment_npoints

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Plot dft and pdos.
        if path_segment_npoints:
            ax.plot(self.dft_eigs[:, 0], label='DFT', color='blue')
            ax.plot(self.dft_eigs, color='blue')
            ax.yaxis.grid(False)  
            ax.set_xticks(
                ticks=np.arange(len(path_special_points))*path_segment_npoints,
                labels=path_special_points,
            )

            # Scatter kpdos on top. 
            for annotated_file_item in annotated_file_list:
                pdos = self.get_kpdos_data(annotated_file_item['filename'])
                ax.scatter(
                    self.dft_eigs,
                    s=pdos,
                    label='{atom} {orbital}'.format(atom=annotated_file_item['atom'], orbital=annotated_file_item['orbital'])
                )
        else:
            xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()    
            ax.plot(xaxis, self.dft_eigs[:, 0], label='DFT', color='blue')
            ax.plot(xaxis, self.dft_eigs, color='blue')
            ax.yaxis.grid(False) 
            ax.set_xticks(
                ticks=special_points,
                labels=special_labels,
            )

             # Scatter kpdos on top. 
            for annotated_file_item in annotated_file_list:
                pdos = self.get_kpdos_data(annotated_file_item['filename'])
                ax.scatter(
                    xaxis,
                    self.dft_eigs,
                    s=pdos,
                    label='{atom} {orbital}'.format(atom=annotated_file_item['atom'], orbital=annotated_file_item['orbital'])
                )

        

        ax.legend()
        ax.set_ylabel('Energy (eV)')
        ax.set_title('k-resolved PDOS')
        if ylim  is not None: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.savefig(save_filename)
        if show: plt.show()
#endregion
