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
class PhbandsPlot:
    def __init__(
        self,
        phbands_filename,
        bandpathpkl_filename,
    ):
        self.phbands_filename = phbands_filename
        self.bandpathpkl_filename = bandpathpkl_filename

        self.num_bands: int = None 
        self.phbands: np.ndarray = None 
        self.kpath: KPath = None 

    def get_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[0]
        self.kpath = load_obj(self.bandpathpkl_filename)
        
    def save_plot(self, save_filename, show=False):
        self.get_data()
        # xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # ax.plot(xaxis, self.phbands, color='blue')
        ax.plot(self.phbands, color='blue')
        ax.yaxis.grid(False)
        # ax.set_xticks(ticks=special_points, labels=special_labels)
        ax.legend()

        fig.savefig(save_filename)

        if show: plt.show()

 #endregion
