#region: Modules.
import matplotlib.pyplot as plt 
import numpy as np 
from fp.io.pkl import load_obj
from fp.structure.kpath import KPath
from fp.flows.fullgridflow import FullGridFlow
from fp.inputs.input_main import Input
import yaml 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PhbandsPlot:
    def __init__(
        self,
        phbands_filename: str = 'struct.freq.gp',
        bandpathpkl_filename: str = 'bandpath.pkl',
    ):
        self.phbands_filename = phbands_filename
        self.bandpathpkl_filename = bandpathpkl_filename

        self.num_bands: int = None 
        self.phbands: np.ndarray = None 
        self.kpath: KPath = None 

    def get_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[1]
        self.kpath = load_obj(self.bandpathpkl_filename)
        
    def save_plot(self, save_filename='phbands.png', show=False, ylim=None):
        # Get some data. 
        self.get_data()
        path_special_points = self.kpath.path_special_points
        path_segment_npoints = self.kpath.path_segment_npoints

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments.
        ax.plot(self.phbands, color='blue')
        ax.yaxis.grid(False)  
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )

        # Set some labels. 
        ax.set_title('Phonon Bandstructure')
        ax.set_ylabel('Freq (cm-1)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.savefig(save_filename)
        if show: plt.show()

class PhonopyPlot(PhbandsPlot):
    def __init__(self, **kwargs):
        super().__init__(phbands_filename='band.yaml', **kwargs)

    def get_data(self):
        with open(self.phbands_filename) as f: data = yaml.safe_load(f)

        nk = len(data['phonon'])
        nb = len(data['phonon'][0]['band'])

        # fill phbands
        self.phbands = np.zeros(shape=(nk, nb), dtype='f8')
        for (k, b), value in np.ndenumerate(self.phbands):
            self.phbands[k, b] = data['phonon'][k]['band'][b]['frequency']*33.356        # Factor in cm^{-1}
 
        self.num_bands = self.phbands.shape[1]
        self.kpath = load_obj(self.bandpathpkl_filename)
 #endregion
