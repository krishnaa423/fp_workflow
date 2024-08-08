#region: Modules.
import numpy as np
import matplotlib.pyplot as plt 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class BseSpectrumPlot:
    def __init__(
        self,
        eh_filename,
        noeh_filename,
    ):
        self.eh_filename = eh_filename
        self.noeh_filename = noeh_filename

    def save_plot(self, save_filename, show=False):
        abs_eh_data = np.loadtxt(self.eh_filename, dtype='f8', skiprows=4)
        abs_noeh_data = np.loadtxt(self.noeh_filename, dtype='f8', skiprows=4)
        
        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(abs_eh_data[:, 0], abs_eh_data[:, 1], label='eh')
        ax.plot(abs_noeh_data[:, 0], abs_noeh_data[:, 1], label='noeh')
        ax.legend()

        fig.savefig(save_filename)

        if show: plt.show()

#endregion
