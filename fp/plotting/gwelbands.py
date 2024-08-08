#region: Modules.
import matplotlib.pyplot as plt 
import numpy as np 
from fp.io.pkl import load_obj
from fp.structure import KPath
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class GwelbandsPlot:
    def __init__(
        self,
        inteqp_filename,
        bandpathpkl_filename,
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
        


    def save_plot(self, save_filename, show=False):
        self.get_data()
        # xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # ax.plot(xaxis, self.emf[:, 0], label='DFT', color='blue')
        # ax.plot(xaxis, self.eqp[:, 0], label='GW', color='green')
        # ax.plot(xaxis, self.emf, color='blue')
        # ax.plot(xaxis, self.eqp, color='green')
        # ax.plot(xaxis, self.emf[:, 0], label='DFT', color='blue')
        # ax.plot(xaxis, self.eqp[:, 0], label='GW', color='green')
        ax.plot(self.emf, color='blue')
        ax.plot(self.eqp, color='green')
        ax.yaxis.grid(False)
        # ax.set_xticks(ticks=special_points, labels=special_labels)
        ax.legend()

        fig.savefig(save_filename)

        if show: plt.show()

 #endregion
