#region: Modules.
from fp.flows import *
from fp.flows.fullgridflow import *
from fp.inputs import *
from fp.inputs.relax import RelaxType
from fp.schedulers import *
from fp.calcs import *
from fp.calcs.relax import Relax
from fp.calcs.wfngeneral import Wfn, Wfnq
from fp.structure import *
from ase import Atoms 
from ase.io import write
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    # Inputs.
    uc_atoms = Atoms(
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
    write('uc_atoms.xsf', uc_atoms)

    fullgridflow = FullGridFlow.from_yml('../fp/data/fullgridflow_sample.yml')
    flowmanage = fullgridflow.get_flowmanage([
            Relax,
            Scf,
            Dfpt,
            Phbands,
            Phdos,
            Phmodes,
            Dos,
            Pdos,
            Dftelbands,
            Kpdos,
            Wannier,
            Wfn,
            Epw,
            Wfnq,
            # Wfnfi,
            # Wfnqfi,
            Epsilon,
            Sigma,
            Inteqp,
            Kernel,
            Absorption,
            Plotxct,
            Bseq,
            # Esf(input=input),
            # Esd(input=input),
            # XctPh(input=input),
            # Pol(input=input),
            # Xctpol(input=input),
        ], 
        save_pkl=True
    )
    flowmanage.create()
    flowmanage.create_job_all_script(
        filename='job_all.sh',
        start_job='job_relax.sh',
        stop_job='job_bseq.sh',
        save_folder_flag=True,
        save_folder='./test_save_folder/',
        flowfile_to_read='flowmanage.pkl'
    )
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion