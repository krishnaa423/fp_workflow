#region: Modules.
from fp.flows import *
from fp.flows.fullgridflow import *
from fp.inputs import *
from fp.inputs.relax import RelaxType
from fp.schedulers import *
from fp.calcs import *
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

    fullgridflow = FullGridFlow(
        scheduler = WSL(),
        
        single_task_desc = JobProcDesc(
            nodes=1, 
            ntasks=1,
            time='00:45', 
        ),
        single_node_desc = JobProcDesc(
            nodes=1,  
            ntasks=1,
            time='00:45',
        ),
        para_desc = JobProcDesc(
            nodes=1, 
            ntasks=1,
            time='00:45', 
        ),
        big_para_desc = JobProcDesc(
            nodes=1, 
            ntasks=1,
            time='00:45', 
        ),
        para_k_desc = JobProcDesc(
            nodes=1,  
            ntasks=1,
            time='00:45',
            # nk=,
        ),
        big_para_k_desc = JobProcDesc(
            nodes=1,  
            ntasks=1,
            time='00:45',
            # nk=,
        ),
        para_epwk_desc = JobProcDesc(
            nodes=1,
            # Tasks should be the same as nk for epwk. 
            ntasks=1,
            time='00:45', 
            # nk=,
        ),
        
        atoms = uc_atoms,
        sc_grid = (1, 1, 1),
        use_esd_atoms_if_needed = True,
        
        path_string = 'LGXWLKG',
        path_npoints = 100,
        
        relax_type = RelaxType.GS_RELAX,

        scf_kgrid=(2, 2, 2),
        scf_cutoff=20.0,

        dfpt_qgrid=(2, 2, 2),
        dfpt_conv_threshold='1.0d-16',
        dfpt_phmode=1,

        dos_kdim = (2, 2, 2) ,
        dftelbands_cond=10,
        wannier_kdim = (2, 2, 2),
        wannier_bands_cond = 10,
        
        wfn_qe_cond = 10,
        wfn_qe_kdim = (2, 2, 2),
        wfn_qe_sym = False,
        wfn_para_cond = 201,
        
        qshift = (0.0, 0.0, 0.001),
        wfnq_qe_kdim = (2, 2, 2),
        wfnq_qe_cond = 10,
        wfnq_qe_sym = False, 
        
        epssig_bands_cond = 200,
        epssig_cutoff = 10.0,
        
        sig_band_val = 4,
        sig_band_cond = 10,
        
        inteqp_band_val = 4,
        
        abs_val_bands = 4,
        abs_cond_bands = 10,
        abs_nevec = 10,

        bseq_Qdim = (2, 2, 2),
        
        plotxct_hole = (0.25, 0.25, 0.25),
        plotxct_sc = (2, 2, 2),
        plotxct_state = 1,
    )

    flow = fullgridflow.get_flowmanage([
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
    ])

    flow.create()
    flow.create_job_all_script(
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