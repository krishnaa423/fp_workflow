#region: Modules.
import numpy as np 
import matplotlib.pyplot as plt 
#endregion

#region: Variables.
#endregion

#region: Functions.
def get_abs():
    eh = np.loadtxt('../absorption_eh.dat', skiprows=4)[:, [0, 1]]
    noeh = np.loadtxt('../absorption_noeh.dat', skiprows=4)[:, [0, 1]]
    return eh, noeh 

def main():
    eh, noeh = get_abs()
    
    print(noeh.shape)
    
    plt.style.use('bmh')
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.plot(eh[:, 0], eh[:, 1], label='eh')
    ax.plot(noeh[:, 0], noeh[:, 1], label='noeh')
    ax.legend()
    ax.set_xlabel('$\omega (eV)$')
    ax.set_ylabel('$\epsilon_2 (Arb)$')
    ax.set_title('Absorption Spectra')
    plt.show()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion