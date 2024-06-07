#region: Modules.
import numpy as np 
import xml.etree.ElementTree as ET 
import matplotlib.pyplot as plt 
import h5py 
#endregion

#region: Variables.
#endregion

#region: Functions.
def gauss(x, mu, sigma):
    # return np.exp(-((x - mu)/sigma)**2)/sigma/np.sqrt(np.pi)
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu)**2 / (2 * sigma**2))

def get_energies():
    root = ET.parse('../dft_elbands.xml').getroot()
    
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
    
    return energies

def get_grids():
    kpath = np.loadtxt('../uc_kpath.txt', dtype='f8')  
    Kpath = np.loadtxt('../sc_Kpath.txt', dtype='f8')  
    Gshift = np.loadtxt('../sc_Gshift.txt', dtype='f8')  
    sc_grid = np.loadtxt('../sc_grid.txt', dtype='f8').astype('i4')
    
    return kpath, Kpath, Gshift, sc_grid

def get_coeff():
    # Set sizes. 
    with h5py.File('../WFN_dft_elbands.h5', 'r') as f:
        fft_grid = f['mf_header/gspace/FFTgrid'][:]
        
        ngk_raw = f['mf_header/kpoints/ngk'][:]
        ngk = np.zeros(shape=(ngk_raw.size+1,))
        ngk[1:] = np.cumsum(ngk_raw)
        
        nrk = f['mf_header/kpoints/nrk'][()] 
        
        coeff_linear = np.vectorize(complex)(f['wfns/coeffs'][:, 0, :, 0], f['wfns/coeffs'][:, 0, :, 1]) # n, kG. 
        
        gvecs = f['wfns/gvecs'][:] # nG, 3. 
        
    
    # Extract data. 
    # print(fft_grid.shape)
    # print(nrk.shape)
    # print(coeff_linear.shape)
    # print(gvecs.shape)
    # print(ngk.shape)
    # print(ngk)
    
    # Fold the grid. 
    gvecs[:, 0] %= fft_grid[0]
    gvecs[:, 1] %= fft_grid[1]
    gvecs[:, 2] %= fft_grid[2]
    
    # c[n, k, Gx, Gy, Gz]. 
    coeff = np.zeros(
        shape=(coeff_linear.shape[0], nrk, fft_grid[0], fft_grid[1], fft_grid[2]),
        dtype='c16',
    )
    
    for band_idx in range(coeff.shape[0]):
        for kpt_idx in range(coeff.shape[1]):
            coeff[
                band_idx, 
                kpt_idx, 
                gvecs[int(ngk[kpt_idx]):int(ngk[kpt_idx+1]), 0], 
                gvecs[int(ngk[kpt_idx]):int(ngk[kpt_idx+1]), 1], 
                gvecs[int(ngk[kpt_idx]):int(ngk[kpt_idx+1]), 2], 
                ] = coeff_linear[band_idx, int(ngk[kpt_idx]):int(ngk[kpt_idx+1])]
    return coeff 

def get_weights(coeff: np.ndarray, Gshift: np.ndarray, sc_grid: np.ndarray):
    
    # coeff: c[n, k, Gx, Gy, Gz]. 
    # weights: w[k, n]. 
    
    weights = np.zeros(shape=(coeff.shape[1], coeff.shape[0]), dtype='f8') 
    
    for kpt_idx in range(weights.shape[0]):
        for band_idx in range(weights.shape[1]):
            weights[kpt_idx, band_idx] = np.sum(np.abs(coeff[band_idx, kpt_idx, 
                                                      int(Gshift[kpt_idx, 0])::int(sc_grid[0]), 
                                                      int(Gshift[kpt_idx, 1])::int(sc_grid[1]), 
                                                      int(Gshift[kpt_idx, 2])::int(sc_grid[2])].flatten()))
            
    return weights

def get_axis_grids():
    elims = np.array([-5, 15])
    ne = 100
    klims = np.array([0, 99])
    nk = 100
    egrid = np.linspace(elims[0], elims[1], ne)
    kgrid = np.linspace(klims[0], klims[1], nk)
    
    return kgrid, egrid 

def get_spectral_map(kgrid, egrid, energies: np.ndarray, weights: np.ndarray):
    sigma = 0.1  # Broadening parameter. 
    
    K, E = np.meshgrid(kgrid, egrid, indexing='ij')
    A = np.zeros_like(K)
    for kpt_idx in range(energies.shape[0]):
        for eig_idx in range(energies.shape[1]):
            mu = energies[kpt_idx, eig_idx]
            A[kpt_idx, :] += weights[kpt_idx, eig_idx]*gauss(egrid, mu, sigma)
            
        # A[kpt_idx, :] /= np.linalg.norm(A[kpt_idx, :])
        
    return K, E, A 

def plot_map(K, E, spectral_map):
    fig = plt.figure()
    ax = fig.add_subplot()
    mesh = ax.pcolormesh(K.transpose(), E.transpose(), spectral_map.transpose(), cmap='magma')
    ax.set_xlabel('Kpoints')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('Spectral function')
    
    fig.colorbar(mesh)

def main():
    
    # Read from dft_elbands.xml -> dft_e. 
    energies = get_energies()
    
    # Read uc_kpath.txt, sc_kpath.txt, sc_gshift.txt as kpath, Kpath, Gshift. 
    kpath, Kpath, Gshift, sc_grid = get_grids()
    
    # Construct c[K, n, Gx, Gy, Gz]. 
    coeff = get_coeff()
    
    # # Construct w[K, n]. 
    weights = get_weights(coeff, Gshift, sc_grid)
    
    # # Construct kgrid, Egrid. 
    kgrid, egrid = get_axis_grids()
    
    # # Construct A[k, E] for kgrid, Egrid.
    K, E, spectral_map = get_spectral_map(kgrid, egrid, energies, weights)
    
    # # Plot A[k, E] and/or save the matrix. 
    plt.style.use('bmh') 
    plot_map(K, E, spectral_map)
    plt.show()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion