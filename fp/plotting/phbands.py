#region: Modules.
import matplotlib.pyplot as plt 
import numpy as np 
from fp.io.pkl import load_obj
from fp.structure import KPath
from fp.flows.fullgridflow import FullGridFlow
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
        fullgridflow_filename,
    ):
        self.phbands_filename = phbands_filename
        self.bandpathpkl_filename = bandpathpkl_filename
        self.fullgridflow_filename = fullgridflow_filename

        self.num_bands: int = None 
        self.phbands: np.ndarray = None 
        self.kpath: KPath = None 
        self.fullgridflow: FullGridFlow = None 

    def get_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[0]
        self.kpath = load_obj(self.bandpathpkl_filename)
        self.fullgridflow = load_obj(self.fullgridflow_filename)
        
    def save_plot(self, save_filename, show=False):
        self.get_data()
        # xaxis, special_points, special_labels = self.kpath.bandpath.get_linear_kpoint_axis()

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()
        path_special_points = self.fullgridflow.path_special_points
        path_segment_npoints = self.fullgridflow.path_segment_npoints

        ax.plot(self.phbands, color='blue')
        ax.yaxis.grid(False)
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )
        sizes = np.zeros_like(self.phbands)
        xaxis = np.arange(sizes.shape[0]).reshape(-1, 1)
        xaxis = np.repeat(xaxis, repeats=sizes.shape[1], axis=1)
        sizes[:, [20, 21, 22, 22, 23, 24]] = np.random.random(size=(sizes.shape[0], 6))*100
        ax.scatter(xaxis, self.phbands, s=sizes, c='green')
        ax.legend()

        fig.savefig(save_filename)

        if show: plt.show()

 #endregion
