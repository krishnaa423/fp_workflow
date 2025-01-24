#region: Modules.
import matplotlib.pyplot as plt 
import numpy as np 
from fp.io.pkl import load_obj
from fp.structure.kpath import KPath
from fp.flows.fullgridflow import FullGridFlow
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class GwelbandsPlot:
    def __init__(
        self,
        inteqp_filename: str ='bandstructure_inteqp.dat',
        bandpathpkl_filename: str ='bandpath.pkl',
    ):
        self.inteqp_filename = inteqp_filename
        self.bandpathpkl_filename = bandpathpkl_filename

        self.num_bands: int = None 
        self.emf: np.ndarray = None 
        self.eqp: np.ndarray = None 
        self.kpath: KPath = None 

    def get_data(self):
        data = np.loadtxt(self.inteqp_filename, skiprows=2)
        num_bands = np.unique(data[:, 1]).size
        emf = data[:, 5].reshape(num_bands, -1).T
        eqp = data[:, 6].reshape(num_bands, -1).T

        self.num_bands = num_bands
        self.emf = emf 
        self.eqp = eqp
        self.kpath = load_obj(self.bandpathpkl_filename)

    def save_plot(self, save_filename='gwelbands.png', show=False, ylim=None, offset=None):
        # Get some data. 
        self.get_data()
        if offset  is not None: self.eqp += offset 
        path_special_points = self.kpath.path_special_points
        path_segment_npoints = self.kpath.path_segment_npoints

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments.
        ax.plot(self.emf[:, 0], label='DFT', color='blue')
        ax.plot(self.eqp[:, 0], label='GW', color='green')
        ax.plot(self.emf, color='blue')
        ax.plot(self.eqp, color='green')
        ax.yaxis.grid(False)  
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )

        ax.set_title('GW Bandstructure')
        ax.set_ylabel('Energy (eV)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        ax.legend()
        fig.savefig(save_filename)
        if show: plt.show()
 #endregion
