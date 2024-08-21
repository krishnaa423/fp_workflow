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
        dryrun_scheduler: Scheduler=None,
        
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
        skip_pseudo_generation: bool = None,
        
        path_special_points: str=None,
        path_segment_npoints: int=None,
        
        relax_type = None,
        read_relaxed_coord: bool = None,

        scf_kgrid=None,
        scf_cutoff=None,

        dfpt_qgrid=None,
        dfpt_conv_threshold:str=None,
        dfpt_phmode: int=None,
        dfpt_njobs: int = None,

        dos_kdim = None ,
        dftelbands_cond: int = None,
        wannier_kdim = None,
        wannier_bands_cond: int = None,
        
        wfn_qe_cond: int = None,
        wfn_qe_kdim = None ,
        wfn_qe_sym = None,
        wfn_para_cond = None,

        epw_exec_loc: str = None,
        
        qshift = None,
        wfnq_qe_cond: int  = None,
        wfnq_qe_kdim = None,
        wfnq_qe_sym = None, 

        wfnfi_qe_cond: int  = None,
        wfnfi_qe_kdim = None,
        wfnfi_qe_sym = None, 

        wfnqfi_qe_cond: int  = None,
        wfnqfi_qe_kdim = None,
        wfnqfi_qe_sym = None, 
        
        epssig_bands_cond: int = None,
        epssig_cutoff: float  = None,
        epssig_wfnlink: str = None,
        epssig_wfnqlink: str = None,
        
        sig_band_val: int = None,
        sig_band_cond: int = None,
        
        inteqp_band_val: int = None,
        inteqp_wfn_co_link: str = None,
        inteqp_wfn_fi_link: str = None,
        
        abs_val_bands: int = None,
        abs_cond_bands: int = None,
        abs_nevec: int = None,
        abs_wfn_co_link: str = None,
        abs_wfnq_co_link: str = None,
        abs_wfn_fi_link: str = None,
        abs_wfnq_fi_link: str = None,

        bseq_Qdim = None,
        
        plotxct_hole = None,
        plotxct_sc = None,
        plotxct_state: int = None,

        esd_fmax: float = None,
        esd_max_steps: int = None,

        xctpol_max_error: float = None,
        xctpol_max_steps: int = None,
    ):
        '''
        Simplifies flow manage. 
        '''
        self.scheduler: Scheduler = scheduler
        self.dryrun_scheduler: Scheduler = dryrun_scheduler
        
        self.single_task_desc: JobProcDesc = JobProcDesc(**single_task_desc)
        self.single_node_desc: JobProcDesc = JobProcDesc(**single_node_desc)
        self.para_desc: JobProcDesc = JobProcDesc(**para_desc)
        self.big_para_desc: JobProcDesc = JobProcDesc(**big_para_desc)
        self.para_k_desc: JobProcDesc = JobProcDesc(**para_k_desc)
        self.big_para_k_desc: JobProcDesc = JobProcDesc(**big_para_k_desc)
        self.para_epwk_desc: JobProcDesc = JobProcDesc(**para_epwk_desc)
        
        self.atoms: str = atoms
        self.sc_grid: np.ndarray = np.array(sc_grid)
        self.use_esd_atoms_if_needed: bool = use_esd_atoms_if_needed
        self.skip_pseudo_generation: bool = skip_pseudo_generation
        
        self.path_special_points: list = path_special_points
        self.path_segment_npoints: int = path_segment_npoints
        
        self.relax_type = relax_type
        self.read_relaxed_coord = read_relaxed_coord

        self.scf_kgrid = scf_kgrid
        self.scf_cutoff = scf_cutoff

        self.dfpt_qgrid = dfpt_qgrid
        self.dfpt_conv_threshold:str = dfpt_conv_threshold
        self.dfpt_phmode: int = dfpt_phmode
        self.dfpt_njobs: int = dfpt_njobs

        self.dos_kdim = dos_kdim 
        self.dftelbands_cond: int = dftelbands_cond
        self.wannier_kdim = wannier_kdim
        self.wannier_bands_cond: int = wannier_bands_cond
        
        self.wfn_qe_cond: int = wfn_qe_cond
        self.wfn_qe_kdim = wfn_qe_kdim 
        self.wfn_qe_sym = wfn_qe_sym
        self.wfn_para_cond = wfn_para_cond

        self.epw_exec_loc: str = epw_exec_loc
        
        self.qshift = qshift
        self.wfnq_qe_cond: int  = wfnq_qe_cond
        self.wfnq_qe_kdim = wfnq_qe_kdim
        self.wfnq_qe_sym = wfnq_qe_sym

        self.wfnfi_qe_cond: int  = wfnfi_qe_cond
        self.wfnfi_qe_kdim = wfnfi_qe_kdim
        self.wfnfi_qe_sym = wfnfi_qe_sym

        self.wfnqfi_qe_cond: int  = wfnqfi_qe_cond
        self.wfnqfi_qe_kdim = wfnqfi_qe_kdim
        self.wfnqfi_qe_sym = wfnqfi_qe_sym
        
        self.epssig_bands_cond: int = epssig_bands_cond
        self.epssig_cutoff: float  = epssig_cutoff
        self.epssig_wfnlink: str = epssig_wfnlink
        self.epssig_wfnqlink: str = epssig_wfnqlink
    
        self.sig_band_val: int = sig_band_val
        self.sig_band_cond: int = sig_band_cond
        
        self.inteqp_band_val: int = inteqp_band_val
        self.inteqp_wfn_co_link: str = inteqp_wfn_co_link
        self.inteqp_wfn_fi_link: str = inteqp_wfn_fi_link
        
        self.abs_val_bands: int = abs_val_bands
        self.abs_cond_bands: int = abs_cond_bands
        self.abs_nevec: int = abs_nevec
        self.abs_wfn_co_link: str = abs_wfn_co_link
        self.abs_wfnq_co_link: str = abs_wfnq_co_link
        self.abs_wfn_fi_link: str = abs_wfn_fi_link
        self.abs_wfnq_fi_link: str = abs_wfnq_fi_link

        self.bseq_Qdim = bseq_Qdim
        
        self.plotxct_hole = plotxct_hole
        self.plotxct_sc = plotxct_sc
        self.plotxct_state: int = plotxct_state

        self.esd_fmax: float = esd_fmax
        self.esd_max_steps: int = esd_max_steps

        self.xctpol_max_error: float = xctpol_max_error
        self.xctpol_max_steps: int = xctpol_max_steps

        # During run. 
        self.max_val: int = None 
        self.input: Input = None 
    
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
            assert hasattr(fullgridflow, key), f'FullGridFlow class does not have attribute: {key}.'

            if 'scheduler' in key:        # Create the scheduler class. 
                first_key, first_value = next(iter(value.items()))
                sched_cls = getattr(schedulers, first_key)
                setattr(fullgridflow, key, sched_cls(**first_value))
            elif '_desc' in key:        # Ones with job descriptions. 
                job_desc = JobProcDesc(**value)
                setattr(fullgridflow, key, job_desc)
            else:
                setattr(fullgridflow, key, value)

        return fullgridflow
                 
    def set_relaxed_coords_from_files(self):
        cell_file = 'relaxed_cell_parameters.txt'
        pos_file = 'relaxed_atomic_positions.txt'

        # Read cell/positions and set them. 
        # Only set if the read files are non-zero in length
        if len(open(cell_file).read())!=0: self.uc_atoms.cell = np.loadtxt(cell_file)
        if len(open(pos_file).read())!=0: self.uc_atoms.positions = np.loadtxt(pos_file)

    def create_atoms(self):
        # Make atoms. 
        self.uc_atoms = read(self.atoms) 

        if self.read_relaxed_coord: self.set_relaxed_coords_from_files()

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
        dryrun = Dryrun(atoms=self.atoms_input, scheduler=self.dryrun_scheduler, job_desc=self.single_node_desc)
        dryrun.create()
        dryrun.run(0.0)
        self.max_val = dryrun.get_max_val()
        dryrun.remove()

    def create_kpath(self):
        self.kpath_obj = KPath(
            atoms=self.uc_atoms,
            path_special_points=self.path_special_points,
            path_segment_npoints=self.path_segment_npoints,
        )
        save_obj(self.kpath_obj, 'bandpath.pkl')
        # self.Kpath, self.Gpath = self.kpath_obj.get_sc_path(self.sc_grid)

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
            njobs=self.dfpt_njobs,
        )

        self.phbands = PhbandsInput(
            kpath=self.kpath_obj,
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
            kpath=self.kpath_obj,
            nbands=self.dftelbands_cond + self.max_val,
            job_desc=self.para_desc,
            job_pw2bgw_desc=self.single_node_desc,
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
            is_reduced=self.wfn_qe_sym,
            bands=self.wfn_qe_cond + self.max_val,
            job_wfn_desc=self.para_k_desc,
            job_pw2bgw_desc=self.single_node_desc,
            job_parabands_desc=self.big_para_desc,
            parabands_bands=self.wfn_para_cond + self.max_val,
        )

        skipped_bands = []
        if self.abs_val_bands!= self.max_val:
            temp = (1, self.max_val - self.abs_val_bands)
            skipped_bands.append(temp)

        if self.abs_cond_bands!= self.wfn_qe_cond:
            temp = (self.max_val + self.abs_cond_bands + 1, self.wfn_qe_cond + self.max_val)
            skipped_bands.append(temp)

        if len(skipped_bands)==0:
            skipped_bands = None

        self.epw = EpwInput(
            kgrid_coarse=self.wfn_qe_kdim,
            qgrid_coarse=self.wfn_qe_kdim,
            kgrid_fine=self.wfn_qe_kdim,
            qgrid_fine=self.wfn_qe_kdim,
            bands=self.abs_cond_bands + self.abs_val_bands,
            exec_loc=self.epw_exec_loc,
            # exec_loc='$SCRATCH/q-e-cpu/bin/epw.x',
            job_desc=self.para_epwk_desc,
            skipped_bands=skipped_bands,     # The input bands are 1 to 14, which are fine.
        )

        self.wfnq = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnq_qe_kdim,
            qshift=self.qshift,
            is_reduced=self.wfnq_qe_sym,
            bands=self.wfnq_qe_cond + self.max_val,
            job_wfn_desc=self.para_k_desc,
            job_pw2bgw_desc=self.single_node_desc,
        )

        self.wfnfi = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnfi_qe_kdim,
            qshift=(0.0, 0.0, 0.000),
            is_reduced=self.wfnfi_qe_sym,
            bands=self.wfnfi_qe_cond,
            job_wfn_desc=self.para_k_desc,
            job_pw2bgw_desc=self.single_node_desc,
        )

        self.wfnqfi = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnqfi_qe_kdim,
            qshift=(0.0, 0.0, 0.001),
            is_reduced=self.wfnqfi_qe_sym,
            bands=self.wfnqfi_qe_cond,
            job_wfn_desc=self.para_k_desc,
            job_pw2bgw_desc=self.single_node_desc,
        )

        self.epsilon = EpsilonInput(
            bands=self.epssig_bands_cond + self.max_val,
            cutoff=self.epssig_cutoff,
            wfn_link=self.epssig_wfnlink,
            wfnq_link=self.epssig_wfnqlink,
            job_desc=self.para_desc,
        )

        self.sigma = SigmaInput(
            bands=self.epssig_bands_cond + self.max_val,
            band_min=self.max_val - self.sig_band_val + 1,
            band_max=self.max_val + self.sig_band_cond,
            cutoff=self.epssig_cutoff,
            wfn_inner_link=self.epssig_wfnlink,
            job_desc=self.para_desc,
        )

        self.inteqp = InteqpInput(
            val_bands_coarse=self.inteqp_band_val,
            cond_bands_coarse=self.dftelbands_cond-1,
            val_bands_fine=self.inteqp_band_val,
            cond_bands_fine=self.dftelbands_cond-1,
            wfn_co_link=self.inteqp_wfn_co_link,
            wfn_fi_link=self.inteqp_wfn_fi_link,
            job_desc=self.para_desc,
        )

        self.kernel = KernelInput(
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            Qshift=(0.0, 0.0, 0.0),
            wfn_co_link=self.abs_wfn_co_link,
            wfnq_co_link=self.abs_wfnq_co_link,
            job_desc=self.para_desc,
        )

        self.absorption = AbsorptionInput(
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            val_bands_fine=self.abs_val_bands,
            cond_bands_fine=self.abs_cond_bands,
            Qshift=(0.0, 0.0, 0.0),
            wfn_co_link=self.abs_wfn_co_link,
            wfnq_co_link=self.abs_wfnq_co_link,
            wfn_fi_link=self.abs_wfn_fi_link,
            wfnq_fi_link=self.abs_wfnq_fi_link,
            num_evec=self.abs_nevec,
            pol_dir=self.qshift,
            job_desc=self.para_desc,
        )
        
        self.plotxct = PlotxctInput(
            hole_position=self.plotxct_hole,
            supercell_size=self.plotxct_sc,
            state=self.plotxct_state,
            wfn_fi_link=self.abs_wfn_fi_link,
            wfnq_fi_link=self.abs_wfnq_fi_link,
            job_desc=self.para_desc,
        )

        self.bseq = BseqInput(
            atoms=self.atoms_input,
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            val_bands_fine=self.abs_val_bands,
            cond_bands_fine=self.abs_cond_bands,
            Qdim=self.bseq_Qdim,
            wfn_co_link=self.abs_wfn_co_link,
            wfnq_co_link=self.abs_wfnq_co_link,
            wfn_fi_link=self.abs_wfn_fi_link,
            wfnq_fi_link=self.abs_wfnq_fi_link,
            num_evec=self.abs_nevec,
            pol_dir=self.qshift,
            job_desc=self.para_desc,
        )

        self.esf = EsfInput(
            job_desc=self.single_task_desc,
        )

        self.esdstep = EsdStepInput(
            fmax=self.esd_fmax,
            max_steps=self.esd_max_steps,
            job_desc=self.single_task_desc,
        )

        self.xctph = XctPhInput(
            job_desc=self.single_task_desc,
        )

        self.xctpol = XctPolInput(
            max_error=self.xctpol_max_error,
            max_steps=self.xctpol_max_steps,
            job_desc=self.single_task_desc,
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
            esf=self.esf,
            esdstep=self.esdstep,
            xctph=self.xctph,
            xctpol=self.xctpol,
        )
        if save: save_obj(self.input, 'input.pkl')

    def create_input(self, save=True):
        
        self.create_atoms()

        if not self.skip_pseudo_generation: self.create_pseudos()

        self.create_kpath()

        self.create_atoms_input()

        self.create_max_val()

        self.create_calcs_input(save)

    def get_flowmanage(self, list_of_step_classes: list, save_pkl: bool =True) -> FlowManage:
        self.create_input(save_pkl)

        list_of_steps = [step_class(self.input) for step_class in list_of_step_classes]
        self.flowmanage: FlowManage = FlowManage(list_of_steps)
        if save_pkl: save_obj(self.flowmanage, 'flowmanage.pkl'); save_obj(self, 'fullgridflow.pkl')
        return self.flowmanage

#endregion
