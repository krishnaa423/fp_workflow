#region: Modules.
import fp 
from fp.flows import *
from fp.inputs import AtomsInput, ScfInput, Input
from fp.schedulers import *
from fp.calcs import *
from ase import Atoms 
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    # Inputs.
    atoms = Atoms(
        symbols=['Si', 'Si'],
        cell=np.array([
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ])*5.43,
        scaled_positions=np.array([
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25],
        ]),
        pbc=[1, 1, 1],
    )

    atoms_input = AtomsInput(atoms=atoms)
    
    scf = ScfInput(
        kdim=(2, 2, 2),
        ecutwfc=20.0,
        job_desc=JobProcDesc(1, 1, '01:00'),
    )

    input = Input(
        scheduler=WSL(),
        atoms=atoms_input,
        scf=scf,
    )

    
    # Flow.
    flow = FlowManage(
        list_of_steps=[
            Scf(input=input),
        ]
    )
    # flow.create()
    flow.remove()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion