#region: Modules.
from ase import Atoms 
import numpy as np 
import os
import subprocess
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
        self.atoms: Atoms = atoms 
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
        pass 

    def get_kpts(self):
        return self.get_ibz_kpts() if self.is_reduced else self.get_fbz_kpts()
#endregion
