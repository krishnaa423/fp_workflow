#region: Modules.
import numpy as np
import matplotlib.pyplot as plt 
from ase.io.jsonio import read_json
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    data = np.loadtxt('../struct.freq.gp')
    kpts = data[:, 0]
    bands = data[:, 1:]
    
    bandpath = read_json('../bandpath.json')
    xcoords, special_xcoords, labels = bandpath.get_linear_kpoint_axis()
    
    plt.style.use('bmh')
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(xcoords, bands, color='black')
    ax.set_xticks(ticks=special_xcoords, labels=labels)
    ax.set_xlabel('Kpoints')
    ax.set_ylabel('$\omega$ $(cm^{-1})$')
    ax.set_xlim(xcoords[0], xcoords[-1])
    ax.set_title('Phonon band structure')
    
    plt.show()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion