#region: Modules.
import h5py 
import numpy as np 
from fp.io import load_obj
from fp.structure import KPath
import matplotlib.pyplot as plt 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPolPlot:
    def __init__(
        self,
        xctpol_filename,
        phbands_filename,
        kpath_filename,
        # Could add projection onto the exciton bands later. 
    ):
        self.xctpol_filename = xctpol_filename
        self.phbands_filename = phbands_filename
        self.kpath_filename = kpath_filename

        # Later. 
        self.xctpol_ph: np.ndarray = None 
        self.kpath: KPath = None 
        self.phbands: np.ndarray = None 

    def get_xctpol(self):
        with h5py.File(self.xctpol_filename, 'r') as r:
            self.xctpol_ph = r['/xctpol_ph'][:]

    def get_phbands(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[0]
        self.kpath: KPath = load_obj(self.kpath_filename)

    def save_plot(self, save_filename, show=False):
        self.get_xctpol()
        self.get_phbands()
        xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        ax.plot(xaxis, self.phbands, color='black', alpha=0.5)
        ax.yaxis.grid(False)
        ax.set_xticks(ticks=special_points, labels=special_labels)

        # Testing for intensity plot. 
        sizes = np.random.random(size=self.phbands.shape)*100
        ax.scatter(np.repeat(xaxis.reshape(-1, 1), repeats=self.phbands.shape[1], axis=1), self.phbands, s=sizes, c='yellow')

        fig.savefig(save_filename)

        if show: plt.show()
#endregion
