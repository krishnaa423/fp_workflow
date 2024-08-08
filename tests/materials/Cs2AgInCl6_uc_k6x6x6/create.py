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
        symbols=['Cs', 'Cs', 'Ag', 'In', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl'],
        cell=np.array([
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ])*10.48,
        scaled_positions=np.array([
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.75],
            [0.50, 0.50, 0.50],
            [0.00, 0.00, 0.00],
            [0.25, 0.25, 0.75],
            [0.25, 0.75, 0.25],
            [0.75, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75],
        ]),
        pbc=[1, 1, 1],
    ) 
    write('uc_atoms.xsf', uc_atoms)

    fullgridflow = FullGridFlow.from_yml('./fullgridflow.yml')
    flowmanage = fullgridflow.get_flowmanage([
            Relax,
            Scf,
            # Dfpt,
            # Phbands,
            # Phdos,
            # Phmodes,
            # Dos,
            # Pdos,
            Dftelbands,
            # Kpdos,
            # Wannier,
            Wfn,
            # Epw,
            Wfnq,
            Wfnfi,
            Wfnqfi,
            Epsilon,
            Sigma,
            Inteqp,
            Kernel,
            Absorption,
            Plotxct,
            # Bseq,
            # Esf,
            # Esd,
            # XctPh,
            # Pol,
            # XctPol,
        ], 
        save_pkl=True
    )
    flowmanage.create()
    fullgridflow.input.scheduler.create_interactive(JobProcDesc(nodes=4, ntasks=16, time="01:45:00"))

    flowmanage.create_job_all_script(
        filename='job_all.sh',
        start_job='job_scf.sh',
        stop_job='job_plotxct.sh',
        flowfile_to_read='flowmanage.pkl',
    )

#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion