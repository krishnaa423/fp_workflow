#region: Modules.
import numpy as np 
import matplotlib.pyplot as plt 
import re 
import glob 
from typing import List, Dict, Tuple
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class PdosPlot:
    def __init__(
        self,
        pdos_prefix: str,
    ):
        self.pdos_prefix: str = pdos_prefix

        # Data.
        pass 

    def get_annotated_file_list(self) -> List[Dict]:
        files = glob.glob(self.pdos_prefix)
        
        annotated_file_list: List[Dict] = []
        for file in files:
            match = re.search(r'.*?_atm#\d+\((?P<atom>.*?)\)_wfc#\d+\((?P<orbital>.*?)\)', file)
            annotated_file_list.append({
                'filename': file,
                'atom': match.group('atom'),
                'orbital': match.group('orbital'),
            })

        return annotated_file_list

    def get_file_data(self, filename) -> Tuple[np.ndarray, np.ndarray]:
        data = np.loadtxt(filename, skiprows=1) 
        energy = data[:, 0] # in eV.
        pdos = data[:, -1] 

        return (energy, pdos)

    def save_plot(self, save_filename, show=False, ylim=None):
        annotated_file_list = self.get_annotated_file_list()

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        
        for annotated_file_item in annotated_file_list:
            energy, pdos = self.get_file_data(annotated_file_item['filename'])
            ax.plot(
                energy,
                pdos,
                label='{atom} {orbital}'.format(atom=annotated_file_item['atom'], orbital=annotated_file_item['orbital'])
            )

        ax.legend()
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS')
        ax.set_title('DOS')
        if ylim  is not None: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.savefig(save_filename)
        if show: plt.show()
#endregion
