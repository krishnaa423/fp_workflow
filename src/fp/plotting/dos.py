#region: Modules.
import numpy as np 
import matplotlib.pyplot as plt 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class DosPlot:
    def __init__(
        self,
        dos_filename: str,
    ):
        self.dos_filename: str = dos_filename

        # Data.
        self.energy: np.ndarray = None  
        self.dos: np.ndarray = None 

    def get_data(self):
        data = np.loadtxt(self.dos_filename, skiprows=1) 
        self.energy = data[:, 0] # in eV.
        self.dos = data[:, 1] 

    def save_plot(self, save_filename, show=False, ylim=None):
        self.get_data()

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        ax.plot(self.energy, self.dos, color='red')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS')
        ax.set_title('DOS')
        if ylim  is not None: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.savefig(save_filename)
        if show: plt.show()
#endregion
