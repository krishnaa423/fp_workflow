#region: Modules.
# TODO: torch, geometric, e3nn, sklearn.
from ase import Atoms
from ase.build import bulk 
from ase.dft.kpoints import get_special_points, bandpath
from ase.io.jsonio import read_json
import numpy as np 
import xml.etree.ElementTree as ET 
import matplotlib.pyplot as plt 

#endregion

#region: Variables.
#endregion

#region: Functions.
def gauss(x, mu, sigma):
    return (1/(sigma * np.sqrt(2 * np.pi))) * np.exp(-(x - mu)**2 / (2 * sigma**2))

def get_energies():
    root = ET.parse('../dft_elbands.xml').getroot()
    
    fermi_energy = float(root.findall('.//fermi_energy')[0].text)*27.2114      # Energy in eV. 
    
    elements = root.findall('.//eigenvalues')
    
    # Set sizes. 
    nk = len(elements)
    nb = len(np.fromstring(elements[0].text, dtype='f8', sep=' '))
    energies = np.zeros(shape=(nk, nb))
    
    # Populate. 
    for kpt_idx, kpt_element in enumerate(elements):
        kpt_energies = np.fromstring(kpt_element.text, dtype='f8', sep=' ')
        energies[kpt_idx, :] = kpt_energies
    energies *= 27.2114 # conversion from Hartree to eV. 
    energies -= fermi_energy    # Offset by the fermi energy. 
    
    return energies

def get_axis_grids():
    elims = np.array([-5, 15])
    ne = 100
    klims = np.array([0, 99])
    nk = 100
    egrid = np.linspace(elims[0], elims[1], ne)
    kgrid = np.linspace(klims[0], klims[1], nk)
    
    return kgrid, egrid 

def get_spectral_map(kgrid, egrid, energies: np.ndarray):
    bandpath = read_json('../bandpath.json')
    xcoords, special_xcoords, labels = bandpath.get_linear_kpoint_axis()
    
    sigma = 0.1  # Broadening parameter. 
    
    K, E = np.meshgrid(xcoords, egrid, indexing='ij')
    A = np.zeros_like(K)
    for kpt_idx in range(energies.shape[0]):
        for eig_idx in range(energies.shape[1]):
            mu = energies[kpt_idx, eig_idx]
            A[kpt_idx, :] += gauss(egrid, mu, sigma)
    
    return K, E, A 

def plot_lines(energies):
    bandpath = read_json('../bandpath.json')
    xcoords, special_xcoords, labels = bandpath.get_linear_kpoint_axis()
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(xcoords, energies, color='black')
    ax.set_xticks(ticks=special_xcoords,labels=labels)
    ax.set_xlabel('Kpoints')
    ax.set_ylabel('Energy (eV)')
    ax.set_xlim(xcoords[0], xcoords[-1])
    ax.set_title('DFT bands lineplot')

def plot_map(K, E, spectral_map):
    bandpath = read_json('../bandpath.json')
    xcoords, special_xcoords, labels = bandpath.get_linear_kpoint_axis()
    
    fig = plt.figure()
    ax = fig.add_subplot()
    mesh = ax.pcolormesh(K.transpose(), E.transpose(), spectral_map.transpose(), cmap='magma')
    ax.set_xticks(ticks=special_xcoords,labels=labels)
    ax.set_xlabel('Kpoints')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('DFT Spectral function')
    
    fig.colorbar(mesh)

def main():
    
    # Create energies matrix. 
    energies = get_energies()
    
    # Create spectral map. 
    K, E, spectral_map = get_spectral_map(*get_axis_grids(), energies)
    
    # Plots. 
    plt.style.use('bmh') 
    plot_lines(energies)
    # plot_map(K, E, spectral_map)
    plt.show()
    
def debug():
    # # Create energies matrix. 
    # energies = get_energies()
    
    # # Create spectral map. 
    # spectral_map: np.ndarray = get_spectral_map(*get_axis_grids(), energies)
    
    # Test gauss. 
    plt.style.use('bmh')
    x = np.linspace(-3, 3, 100)
    y = gauss(x, 0, 0.1)
    plt.plot(x, y)
    plt.show()
    
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
    # debug()
#endregion