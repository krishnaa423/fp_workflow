#region: Modules.
from ase import Atoms 
import numpy as np 
import os
import subprocess
from fp.inputs import *
from io import StringIO
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Kgrid:
    def __init__(
        self, 
        atoms,
        kdim,
        qshift=(0.0, 0.0, 0.0),
        is_reduced=False,
    ):
        self.kdim: np.ndarray = np.array(kdim).astype(dtype='i4')
        self.atoms: AtomsInput = atoms 
        self.qshift: np.ndarray = np.array(qshift)
        self.is_reduced: bool = is_reduced

    def get_fbz_kpts(self):
        # Calc the kpts.
        command = ['kmesh.pl']
        args = [f'{int(self.kdim[0])}', f'{int(self.kdim[1])}', f'{int(self.kdim[2])}']
        result = subprocess.run(command + args, capture_output=True, text=True)

        text_io = StringIO('\n'.join(result.stdout.splitlines()[2:]))

        # Read the kpts.
        kpts = np.loadtxt(text_io, dtype='f8')

        # Reshape if needed.
        if kpts.ndim == 1 : kpts = kpts.reshape(1, kpts.size)

        # Set the last column to one.
        kpts[:, 3] = 1.0

        # Add the qshift.
        kpts[:, 0] += self.qshift[0]
        kpts[:, 1] += self.qshift[1]
        kpts[:, 2] += self.qshift[2]

        return kpts

    def get_ibz_kpts(self):
        with open('kgrid.inp', 'w') as f:
            f.write(f'{self.kdim[0]} {self.kdim[1]} {int(self.kdim[2])}\n')     
            f.write(f'0.0 0.0 0.0\n')     
            f.write(f'{self.qshift[0]:15.10f} {self.qshift[1]:15.10f} {self.qshift[2]:15.10f}\n')
            f.write(f'{self.atoms.get_scf_cell()}')
            f.write(f'{self.atoms.get_nat()}\n')
            f.write(f'{self.atoms.get_scf_atomic_positions(first_column="atom_index")}')
            f.write(f'0 0 0\n')
            f.write(f'.false.\n')
            f.write(f'.false.\n')
            f.write(f'.false.\n')
        
        command = ['kgrid.x']
        args = ['kgrid.inp', 'kgrid.log', 'kgrid.out']
        result = subprocess.run(command + args, capture_output=True, text=True)
        
        kpts = np.loadtxt('kgrid.log', skiprows=2)
        
        if kpts.ndim == 1 : kpts = kpts.reshape(1, kpts.size)

        return kpts 

    def get_kpts(self):
        return self.get_ibz_kpts() if self.is_reduced else self.get_fbz_kpts()
#endregion
