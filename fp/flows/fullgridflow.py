#region: Modules.
from fp.flows import *
from fp.flows.flow_manage import *
from fp.inputs import *
from fp.schedulers import *
import fp.schedulers as schedulers
from fp.calcs import *
from fp.calcs.dryrun import *
from fp.structure import *
from ase import Atoms 
from ase.io import write, read
import numpy as np 
from ase.build import make_supercell
from fp.io import *
import yaml 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class FullGridFlow:
    def __init__(
        self,
        scheduler: Scheduler=None,
        
        single_task_desc: dict=dict(),
        single_node_desc: dict=dict(),
        para_desc: dict=dict(),
        big_para_desc: dict=dict(),
        para_k_desc: dict=dict(),
        big_para_k_desc: dict=dict(),
        para_epwk_desc: dict=dict(),
        
        atoms: str=None,
        sc_grid: np.ndarray=None,
        use_esd_atoms_if_needed: bool = None,
        
        path_string: str=None,
        path_npoints: int=None,
        
        relax_type = None,

        scf_kgrid=None,
        scf_cutoff=None,

        dfpt_qgrid=None,
        dfpt_conv_threshold:str=None,
        dfpt_phmode: int=None,

        dos_kdim = None ,
        dftelbands_cond: int = None,
        wannier_kdim = None,
        wannier_bands_cond: int = None,
        
        wfn_qe_cond: int = None,
        wfn_qe_kdim = None ,
        wfn_qe_sym = None,
        wfn_para_cond = None,
        
        qshift = None,
        wfnq_qe_cond: int  = None,
        wfnq_qe_kdim = None,
        wfnq_qe_sym = None, 
        
        epssig_bands_cond: int = None,
        epssig_cutoff: float  = None,
        
        sig_band_val: int = None,
        sig_band_cond: int = None,
        
        inteqp_band_val: int = None,
        
        abs_val_bands: int = None,
        abs_cond_bands: int = None,
        abs_nevec: int = None,

        bseq_Qdim = None,
        
        plotxct_hole = None,
        plotxct_sc = None,
        plotxct_state: int = None,
    ):
        '''
        Simplifies flow manage. 
        '''
        self.scheduler: Scheduler=scheduler
        
        self.single_task_desc: JobProcDesc=JobProcDesc(**single_task_desc)
        self.single_node_desc: JobProcDesc=JobProcDesc(**single_node_desc)
        self.para_desc: JobProcDesc=JobProcDesc(**para_desc)
        self.big_para_desc: JobProcDesc=JobProcDesc(**big_para_desc)
        self.para_k_desc: JobProcDesc=JobProcDesc(**para_k_desc)
        self.big_para_k_desc: JobProcDesc=JobProcDesc(**big_para_k_desc)
        self.para_epwk_desc: JobProcDesc=JobProcDesc(**para_epwk_desc)
        
        self.atoms: str=atoms
        self.sc_grid: np.ndarray=np.array(sc_grid)
        self.use_esd_atoms_if_needed: bool = use_esd_atoms_if_needed
        
        self.path_string: str=path_string
        self.path_npoints: int=path_npoints
        
        self.relax_type = relax_type

        self.scf_kgrid=scf_kgrid
        self.scf_cutoff=scf_cutoff

        self.dfpt_qgrid=dfpt_qgrid
        self.dfpt_conv_threshold:str=dfpt_conv_threshold
        self.dfpt_phmode: int=dfpt_phmode

        self.dos_kdim = dos_kdim 
        self.dftelbands_cond: int = dftelbands_cond
        self.wannier_kdim = wannier_kdim
        self.wannier_bands_cond: int = wannier_bands_cond
        
        self.wfn_qe_cond: int = wfn_qe_cond
        self.wfn_qe_kdim = wfn_qe_kdim 
        self.wfn_qe_sym = wfn_qe_sym
        self.wfn_para_cond = wfn_para_cond
        
        self.qshift = qshift
        self.wfnq_qe_cond: int  = wfnq_qe_cond
        self.wfnq_qe_kdim = wfnq_qe_kdim
        self.wfnq_qe_sym = wfnq_qe_sym
        
        self.epssig_bands_cond: int = epssig_bands_cond
        self.epssig_cutoff: float  = epssig_cutoff
    
        self.sig_band_val: int = sig_band_val
        self.sig_band_cond: int = sig_band_cond
        
        self.inteqp_band_val: int = inteqp_band_val
        
        self.abs_val_bands: int = abs_val_bands
        self.abs_cond_bands: int = abs_cond_bands
        self.abs_nevec: int = abs_nevec

        self.bseq_Qdim = bseq_Qdim
        
        self.plotxct_hole = plotxct_hole
        self.plotxct_sc = plotxct_sc
        self.plotxct_state: int = plotxct_state
    
    @staticmethod
    def from_yml(filename):
        '''
        Generate a fullgrid flow object from a yml file.
        '''
        # Open and read the YAML file
        with open(filename, 'r') as file:
            data: dict = yaml.safe_load(file)

        fullgridflow: FullGridFlow = FullGridFlow()
        for key, value in data.items():
            # Debugging.
            # print(f'key: {key}, value: {value}')

            if key=='scheduler':
                sched_cls = getattr(schedulers, value)
                setattr(fullgridflow, key, sched_cls())
            elif key in ['single_task_desc', 'single_node_desc', 'big_para_desc', 'para_k_desc', 'big_para_k_desc', 'para_epwk_desc']:
                setattr(fullgridflow, key, JobProcDesc(value))
            else:
                setattr(fullgridflow, key, value)

        return fullgridflow
         
    def create_atoms(self):
        # Make atoms. 
        self.uc_atoms = read(self.atoms) 
        self.sc_atoms = make_supercell(self.uc_atoms, np.diag(self.sc_grid))

        # Replace with ESD atoms if needed. 
        if self.use_esd_atoms_if_needed:
            if os.path.exists('./esd_atoms.xsf'): 
                self.sc_atoms = read('./esd_atoms.xsf')

        # Save XSF structure files.
        write('uc_atoms.xsf', self.uc_atoms)
        write('sc_atoms.xsf', self.sc_atoms)

    def create_pseudos(self):
        FlowManage.create_pseudos(self.uc_atoms, is_fr=False)

    def create_atoms_input(self):
        self.atoms_input = AtomsInput(atoms=self.sc_atoms)

    def create_max_val(self):
        dryrun = Dryrun(atoms=self.atoms_input)
        dryrun.create()
        dryrun.run(0.0)
        self.max_val = dryrun.get_max_val()
        dryrun.remove()

    def create_kpath(self):
        self.kpath_obj = KPath(
            atoms=self.uc_atoms,
            path_string=self.path_string,
            npoints=100,
        )
        save_obj(self.kpath_obj, 'bandpath.pkl')
        self.Kpath, self.Gpath = self.kpath_obj.get_sc_path(self.sc_grid)

    def create_calcs_input(self, save=True):
        self.relax = RelaxInput(
            max_val=self.max_val,
            job_desc=self.para_desc,
            relax_type=self.relax_type,
        )

        self.scf = ScfInput(
            kdim=self.scf_kgrid,
            ecutwfc=self.scf_cutoff,
            job_desc=self.para_desc,
        )

        self.dfpt = DfptInput(
            atoms=self.atoms_input,
            qgrid=self.dfpt_qgrid,
            conv_threshold=self.dfpt_conv_threshold,
            job_desc=self.para_k_desc,
        )

        self.phbands = PhbandsInput(
            kpath=self.Kpath,
            job_desc=self.para_k_desc,
        )

        self.phdos = PhdosInput(
            qdim=self.dos_kdim,
            job_desc=self.para_k_desc,
        )

        self.phmodes = PhmodesInput(
            qidx=self.dfpt_phmode,
            job_desc=self.para_desc,
        )

        self.dos = DosInput(
            kdim=self.dos_kdim,
            bands=self.dftelbands_cond + self.max_val,
            job_desc=self.para_desc,
        )

        self.dftelbands = DftelbandsInput(
            kpath=self.Kpath,
            nbands=self.dftelbands_cond + self.max_val,
            job_desc=self.para_desc,
            job_pw2bgw_desc=self.para_desc,
        )

        self.kpdos = KpdosInput(
            job_desc = self.para_desc,
        )

        self.wannier = WannierInput(
            atoms=self.atoms_input,
            kdim=self.wannier_kdim,
            num_bands=self.wannier_bands_cond + self.max_val,
            num_wann=self.wannier_bands_cond + self.max_val,
            job_wfnwan_desc=self.para_desc,
            job_pw2wan_desc=self.single_node_desc,
            job_wan_desc=self.para_epwk_desc,
        )

        self.wfn = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfn_qe_kdim,
            qshift=(0.0, 0.0, 0.0),
            is_reduced=False,
            bands=self.wfn_qe_cond + self.max_val,
            job_wfn_desc=self.para_k_desc,
            job_pw2bgw_desc=self.single_node_desc,
            job_parabands_desc=self.big_para_desc,
            parabands_bands=self.wfn_para_cond + self.max_val,
        )

        if self.abs_val_bands==self.max_val:
            if self.wfn_qe_cond==self.abs_cond_bands:
                skipped_bands=None
            else:
                skipped_bands = [(self.max_val + self.abs_cond_bands + 1, self.wfn_qe_cond + self.max_val)]
        else:
            skipped_bands = [(1, self.max_val - self.abs_val_bands), (self.max_val + self.abs_cond_bands + 1, self.wfn_qe_cond + self.max_val)]
        self.epw = EpwInput(
            kgrid_coarse=self.wfn_qe_kdim,
            qgrid_coarse=self.wfn_qe_kdim,
            kgrid_fine=self.wfn_qe_kdim,
            qgrid_fine=self.wfn_qe_kdim,
            bands=self.wfnq_qe_cond + self.max_val,
            exec_loc='$SCRATCH/q-e-cpu/bin/epw.x',
            job_desc=self.para_epwk_desc,
            skipped_bands=skipped_bands,     # The input bands are 1 to 14, which are fine.
        )

        self.wfnq = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnq_qe_kdim,
            qshift=self.qshift,
            is_reduced=False,
            bands=self.wfnq_qe_cond + self.max_val,
            job_wfn_desc=self.para_k_desc,
            job_pw2bgw_desc=self.single_node_desc,
        )

        # self.wfnfi = WfnGeneralInput(
        #     atoms=self.atoms_input,
        #     kdim=(2, 2, 2),
        #     qshift=(0.0, 0.0, 0.000),
        #     is_reduced=False,
        #     bands=14,
        #     job_wfn_desc=self.para_k_desc,
        #     job_pw2bgw_desc=self.single_node_desc,
        # )

        # self.wfnqfi = WfnGeneralInput(
        #     atoms=self.atoms_input,
        #     kdim=(2, 2, 2),
        #     qshift=(0.0, 0.0, 0.001),
        #     is_reduced=False,
        #     bands=14,
        #     job_wfn_desc=self.para_k_desc,
        #     job_pw2bgw_desc=self.single_node_desc,
        # )

        self.epsilon = EpsilonInput(
            bands=self.epssig_bands_cond + self.max_val,
            cutoff=self.epssig_cutoff,
            wfn_link='WFN_parabands.h5',
            wfnq_link='WFNq_coo.h5',
            job_desc=self.para_desc,
        )

        self.sigma = SigmaInput(
            bands=self.epssig_bands_cond + self.max_val,
            band_min=self.max_val - self.sig_band_val + 1,
            band_max=self.max_val + self.sig_band_cond,
            cutoff=self.epssig_cutoff,
            wfn_inner_link='WFN_parabands.h5',
            job_desc=self.para_desc,
        )

        self.inteqp = InteqpInput(
            val_bands_coarse=self.inteqp_band_val,
            cond_bands_coarse=self.dftelbands_cond-1,
            val_bands_fine=self.inteqp_band_val,
            cond_bands_fine=self.dftelbands_cond-1,
            wfn_co_link='./WFN_coo',
            wfn_fi_link='WFN_dftelbands',
            job_desc=self.para_desc,
        )

        self.kernel = KernelInput(
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            Qshift=(0.0, 0.0, 0.0),
            wfn_co_link='./WFN_parabands.h5',
            wfnq_co_link='./WFN_parabands.h5',
            job_desc=self.para_desc,
        )

        self.absorption = AbsorptionInput(
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            val_bands_fine=self.abs_val_bands,
            cond_bands_fine=self.abs_cond_bands,
            Qshift=(0.0, 0.0, 0.0),
            wfn_co_link='./WFN_parabands.h5',
            wfnq_co_link='./WFN_parabands.h5',
            wfn_fi_link='./WFN_parabands.h5',
            wfnq_fi_link='./WFN_parabands.h5',
            num_evec=self.abs_nevec,
            pol_dir=self.qshift,
            job_desc=self.para_desc,
        )
        
        self.plotxct = PlotxctInput(
            hole_position=self.plotxct_hole,
            supercell_size=self.plotxct_sc,
            state=self.plotxct_state,
            wfn_fi_link='./WFN_parabands.h5',
            wfnq_fi_link='./WFNq_coo.h5',
            job_desc=self.para_desc,
        )

        self.bseq = BseqInput(
            atoms=self.atoms_input,
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            val_bands_fine=self.abs_val_bands,
            cond_bands_fine=self.abs_cond_bands,
            Qdim=self.bseq_Qdim,
            wfn_co_link='./WFN_parabands.h5',
            wfnq_co_link='./WFN_parabands.h5',
            wfn_fi_link='./WFN_parabands.h5',
            wfnq_fi_link='./WFN_parabands.h5',
            num_evec=self.abs_nevec,
            pol_dir=self.qshift,
            job_desc=self.para_desc,
        )

        self.input: Input = Input(
            scheduler=self.scheduler,
            atoms=self.atoms_input,
            scf=self.scf,
            relax=self.relax,
            dfpt=self.dfpt,
            phbands=self.phbands,
            phdos=self.phdos,
            phmodes=self.phmodes,
            dos=self.dos,
            dftelbands=self.dftelbands,
            kpdos=self.kpdos,
            wannier=self.wannier,
            wfn=self.wfn,
            epw=self.epw,
            wfnq=self.wfnq,
            wfnfi=self.wfn,
            wfnqfi=self.wfnq,
            epsilon=self.epsilon,
            sigma=self.sigma,
            inteqp=self.inteqp,
            kernel=self.kernel,
            absorption=self.absorption,
            plotxct=self.plotxct,
            bseq=self.bseq,
        )
        if save: save_obj(self.input, 'input.pkl')

    def create_input(self, save=True):
        
        self.create_atoms()

        self.create_pseudos()

        self.create_kpath()

        self.create_atoms_input()

        self.create_max_val()

        self.create_calcs_input(save)

    def get_flowmanage(self, list_of_step_classes: list, save_pkl: bool =True) -> FlowManage:
        self.create_input(save_pkl)

        list_of_steps = [step_class(self.input) for step_class in list_of_step_classes]
        self.flowmanage: FlowManage = FlowManage(list_of_steps)
        if save_pkl: save_obj(self.flowmanage, 'flowmanage.pkl')
        return self.flowmanage

#endregion
