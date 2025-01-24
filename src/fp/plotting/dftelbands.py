#region: Modules.
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt 
import numpy as np 
from fp.io.pkl import load_obj
from fp.structure.kpath import KPath
from fp.flows.fullgridflow import FullGridFlow
from ase.units import Hartree, eV
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DftelbandsPlot:
    def __init__(
        self,
        dftelbands_xml_filename: str ='dftelbands.xml',
        bandpathpkl_filename: str = 'bandpath.pkl',
    ):
        self.dftelbands_xml_filename = dftelbands_xml_filename
        self.bandpathpkl_filename = bandpathpkl_filename

        self.num_bands: int = None 
        self.emf: np.ndarray = None 
        self.kpath: KPath = None 
        self.dft_eigs: np.ndarray = None 
    
    def get_data(self):
        tree = ET.parse(self.dftelbands_xml_filename)
        root = tree.getroot()

        eig_nodes = root.findall('.//ks_energies/eigenvalues')
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree
        num_kpts = len(eig_nodes)
        num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
        dft_eigs = np.zeros(shape=(num_kpts, num_bands), dtype='f8')
        for kpt_idx, node in enumerate(eig_nodes):
            dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*Hartree - fermi_energy

        self.dft_eigs = dft_eigs

        self.kpath = load_obj(self.bandpathpkl_filename)

    def save_plot(self, save_filename='dftelbands.png', show=False, ylim=None):
        self.get_data()
        path_special_points = self.kpath.path_special_points
        path_segment_npoints = self.kpath.path_segment_npoints

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments or total npoints. 
        ax.plot(self.dft_eigs[:, 0], label='DFT', color='blue')
        ax.plot(self.dft_eigs, color='blue')
        ax.yaxis.grid(False)  
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )
        # else:
        #     xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()    
        #     ax.plot(xaxis, self.dft_eigs[:, 0], label='DFT', color='blue')
        #     ax.plot(xaxis, self.dft_eigs, color='blue')
        #     ax.yaxis.grid(False) 
        #     ax.set_xticks(
        #         ticks=special_points,
        #         labels=special_labels,
        #     )

        ax.set_title('DFT Bandstructure')
        ax.set_ylabel('Energy (eV)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        ax.legend()
        fig.savefig(save_filename)
        if show: plt.show()
#endregion