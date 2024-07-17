#region: Modules.
from fp.flows import *
from fp.inputs import *
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

    kpath = atoms.cell.bandpath('LGXWLKG', 100).kpts

    job_desc_node = JobProcDesc(nodes=1, ntasks=1, time='01:00')

    atoms_input = AtomsInput(atoms=atoms)

    FlowManage.create_pseudos(atoms, is_fr=False)
    dryrun = Dryrun(atoms=atoms_input)
    dryrun.create()
    dryrun.run(0.0)
    max_val = dryrun.get_max_val()
    dryrun.remove()
    
    scf = ScfInput(
        kdim=(2, 2, 2),
        ecutwfc=20.0,
        job_desc=job_desc_node,
    )

    relax = RelaxInput(
        max_val=max_val,
        job_desc=job_desc_node,
    )

    dfpt = DfptInput(
        atoms=atoms_input,
        qgrid=(2, 2, 2),
        conv_threshold='1.0d-16',
        job_desc=job_desc_node,
    )

    phbands = PhbandsInput(
        kpath=kpath,
        job_desc=job_desc_node,
    )

    phdos = PhdosInput(
        qdim=(2, 2, 2),
        job_desc=job_desc_node,
    )

    phmodes = PhmodesInput(
        qidx=1,
        job_desc=job_desc_node,
    )

    dos = DosInput(
        kdim=(2, 2, 2),
        bands=14,
        job_desc=job_desc_node,
    )

    dftelbands = DftelbandsInput(
        kpath=kpath,
        nbands=14,
        job_desc=job_desc_node,
        job_pw2bgw_desc=job_desc_node,
    )

    kpdos = KpdosInput(
        job_desc = job_desc_node,
    )

    wannier = WannierInput(
        atoms=atoms_input,
        kdim=(2, 2, 2),
        num_bands=14,
        num_wann=14,
        job_wfnwan_desc=job_desc_node,
        job_pw2wan_desc=job_desc_node,
        job_wan_desc=job_desc_node,
    )

    wfn = WfnGeneralInput(
        atoms=atoms_input,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.0),
        is_reduced=False,
        bands=14,
        job_wfn_desc=job_desc_node,
        job_pw2bgw_desc=job_desc_node,
        job_parabands_desc=job_desc_node,
        parabands_bands=201,
    )

    epw = EpwInput(
        kgrid_coarse=(2, 2, 2),
        qgrid_coarse=(2, 2, 2),
        kgrid_fine=(2, 2, 2),
        qgrid_fine=(2, 2, 2),
        bands=14,
        exec_loc='$SCRATCH/q-e-cpu/bin/epw.x',
        job_desc=job_desc_node,
        skipped_bands=None,     # The input bands are 1 to 14, which are fine.
    )

    wfnq = WfnGeneralInput(
        atoms=atoms_input,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.001),
        is_reduced=False,
        bands=14,
        job_wfn_desc=job_desc_node,
        job_pw2bgw_desc=job_desc_node,
    )

    wfnfi = WfnGeneralInput(
        atoms=atoms_input,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.000),
        is_reduced=False,
        bands=14,
        job_wfn_desc=job_desc_node,
        job_pw2bgw_desc=job_desc_node,
    )

    wfnqfi = WfnGeneralInput(
        atoms=atoms_input,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.001),
        is_reduced=False,
        bands=14,
        job_wfn_desc=job_desc_node,
        job_pw2bgw_desc=job_desc_node,
    )

    epsilon = EpsilonInput(
        bands=200,
        cutoff=10.0,
        wfn_link='WFN_parabands.h5',
        wfnq_link='WFNq_coo.h5',
        job_desc=job_desc_node,
    )

    sigma = SigmaInput(
        bands=200,
        band_min=1,
        band_max=14,
        cutoff=10.0,
        wfn_inner_link='WFN_parabands.h5',
        job_desc=job_desc_node,
    )

    input = Input(
        scheduler=WSL(),
        atoms=atoms_input,
        scf=scf,
        relax=relax,
        dfpt=dfpt,
        phbands=phbands,
        phdos=phdos,
        phmodes=phmodes,
        dos=dos,
        dftelbands=dftelbands,
        kpdos=kpdos,
        wannier=wannier,
        wfn=wfn,
        epw=epw,
        wfnq=wfnq,
        wfnfi=wfnfi,
        wfnqfi=wfnqfi,
        epsilon=epsilon,
        sigma=sigma,
    )

    # Flow.
    flow = FlowManage(
        list_of_steps=[
            # Relax(input=input),
            Scf(input=input),
            # Dfpt(input=input),
            # Phbands(input=input),
            # Phdos(input=input),
            # Phmodes(input=input),
            # Dos(input=input),
            # Pdos(input=input),
            # Dftelbands(input=input),
            # Kpdos(input=input),
            # Wannier(input=input),
            # Wfngeneral(input=input),
            Wfn(input=input),
            # Epw(input=input),
            Wfnq(input=input),
            # Wfnfi(input=input),
            # Wfnqfi(input=input),
            Epsilon(input=input),
            Sigma(input=input),
            # Inteqp(input=input),
            # Kernel(input=input),
            # Absorption(input=input),
            # Plotxct(input=input),
            # BseQ(input=input),
            # XctPh(input=input),
            # Esf(input=input),
            # Esd(input=input),
            # Pol(input=input),
            # Xctpol(input=input),
        ]
    )
    flow.create()
    flow.run()
    # flow.remove()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion