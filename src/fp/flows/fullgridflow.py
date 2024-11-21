#region: Modules.
from fp.flows import *
from fp.flows.flow_manage import *
from fp.inputs import *
from fp.schedulers.scheduler import *
import fp.schedulers as schedulers
from fp.jobs import *
from fp.jobs.dryrun import *
from fp.structure import *
from ase import Atoms 
from ase.io import write, read
import numpy as np 
from ase.build import make_supercell
from fp.io.pkl import *
import yaml 
from typing import List
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class FullGridFlow:
    def __init__(
        self,
        **kwargs,
    ):
        #region: Class attributes. 
        self.scheduler: Scheduler = None
        self.scheduler_dryrun: Scheduler = None
        
        self.job_single_task: JobProcDesc = None
        self.job_single_node: JobProcDesc = None
        self.job_para: JobProcDesc = None
        self.job_big_para: JobProcDesc = None
        self.job_para_k: JobProcDesc = None
        self.job_big_para_k: JobProcDesc = None
        self.job_para_epwk: JobProcDesc = None
        
        self.atoms: str = None
        self.atoms_sc_grid: List[int] = None
        self.atoms_use_esd_if_needed: bool = True 
        self.atoms_skip_pseudo_generation: bool = False 

        self.path_special_points: list = None
        self.path_segment_npoints: int = None
        self.path_total_npoints: int = None
        
        self.relax_type: RelaxType = None 
        self.relax_read_coord: bool = False  
        self.relax_extra_control_args: str = None  
        self.relax_extra_system_args: str = None  
        self.relax_extra_electrons_args: str = None  

        self.scf_kgrid: List[int] = None 
        self.scf_cutoff: float = None
        self.scf_is_spinorbit: bool = False
        self.scf_extra_control_args: str = None  
        self.scf_extra_system_args: str = None  
        self.scf_extra_electrons_args: str = None  

        self.dfpt_qgrid: List[int] = None 
        self.dfpt_conv_threshold:str = None
        self.dfpt_extra_args: str = None 
        self.phbands_extra_q2r_args: str = None 
        self.phbands_extra_matdyn_args: str = None 
        self.phdos_extra_q2r_args: str = None 
        self.phdos_extra_matdyn_args: str = None 
        self.phmodes_extra_args: str = None
        self.phmodes_qpt_idx: int = None # Starts from 1.  

        self.dos_kdim: List[int] = None 
        self.dos_extra_control_args: str = None 
        self.dos_extra_system_args: str = None 
        self.dos_extra_electrons_args: str = None
        self.dos_extra_args: str = None  
        self.pdos_extra_args: str = None   
        self.dftelbands_cond: int = None
        self.dftelbands_extra_control_args: str = None  
        self.dftelbands_extra_system_args: str = None  
        self.dftelbands_extra_electrons_args: str = None
        self.kpdos_extra_args: str = None   
        self.wannier_kdim: List[int] = None
        self.wannier_bands_cond: int = None
        self.wannier_extra_control_args: str = None 
        self.wannier_extra_system_args: str = None 
        self.wannier_extra_electrons_args: str = None 
        self.wannier_extra_args: str = None 
        
        self.wfn_qe_cond: int = None
        self.wfn_qe_kdim: List[int] = None 
        self.wfn_qe_sym: bool = None
        self.wfn_para_cond: int = None
        self.wfn_extra_control_args: str = None 
        self.wfn_extra_system_args: str = None 
        self.wfn_extra_electrons_args: str = None 
        self.wfn_extra_parabands_args: str = None

        self.epw_exec_loc: str = None
        self.epw_extra_args: str = None 
        
        self.qshift: List[float] = None
        self.wfnq_qe_cond: int  = None
        self.wfnq_qe_kdim: List[int] = None
        self.wfnq_qe_sym = None
        self.wfnq_extra_control_args: str = None 
        self.wfnq_extra_system_args: str = None 
        self.wfnq_extra_electrons_args: str = None 

        self.wfnfi_qe_cond: int  = None
        self.wfnfi_qe_kdim: List[int] = None
        self.wfnfi_qe_sym: List[int] = None
        self.wfnfi_extra_control_args: str = None 
        self.wfnfi_extra_system_args: str = None 
        self.wfnfi_extra_electrons_args: str = None 

        self.wfnqfi_qe_cond: int  = None
        self.wfnqfi_qe_kdim: List[int] = None
        self.wfnqfi_qe_sym: bool = None
        self.wfnqfi_extra_control_args: str = None 
        self.wfnqfi_extra_system_args: str = None 
        self.wfnqfi_extra_electrons_args: str = None 
        
        self.epssig_bands_cond: int = None
        self.epssig_cutoff: float  = None
        self.epssig_wfnlink: str = None
        self.epssig_wfnqlink: str = None
        self.eps_extra_args: str = None
    
        self.sig_band_val: int = None
        self.sig_band_cond: int = None
        self.sig_extra_args: str = None
        
        self.inteqp_band_val: int = None
        self.inteqp_wfn_co_link: str = None
        self.inteqp_wfn_fi_link: str = None
        self.inteqp_extra_args: str = None 
        
        self.abs_val_bands: int = None
        self.abs_cond_bands: int = None
        self.abs_nevec: int = None
        self.abs_wfn_co_link: str = None
        self.abs_wfnq_co_link: str = None
        self.abs_wfn_fi_link: str = None
        self.abs_wfnq_fi_link: str = None
        self.abs_pol_dir: list = None
        self.ker_extra_args: str = None
        self.abs_extra_args: str = None

        self.bseq_Qdim: List[int] = None
        
        self.plotxct_hole: List[float] = None
        self.plotxct_sc: List[int] = None
        self.plotxct_state: int = None
        self.plotxct_extra_args: str = None 
        self.plotxct_wfn_fi_link: str = None
        self.plotxct_wfnq_fi_link: str = None 

        self.xctpol_max_error: float = None
        self.xctpol_max_steps: int = None

        # During run. 
        self.max_val: int = None 
        self.input: Input = None 
        #endregion
    
        for key, value in kwargs.items():
            setattr(self, key, value)

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
            elif 'job_' in key:        # Ones with job descriptions. 
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

        if self.relax_read_coord: self.set_relaxed_coords_from_files()

        self.sc_atoms = make_supercell(self.uc_atoms, np.diag(self.atoms_sc_grid))

        # Replace with ESD atoms if needed. 
        if self.atoms_use_esd_if_needed:
            if os.path.exists('./esd_atoms.xsf'): 
                self.sc_atoms = read('./esd_atoms.xsf')

        # Save XSF structure files.
        write('uc_atoms.xsf', self.uc_atoms)
        write('sc_atoms.xsf', self.sc_atoms)

    def create_pseudos(self):
        FlowManage.create_pseudos(self.uc_atoms, is_fr=self.scf_is_spinorbit)

    def create_atoms_input(self):
        self.atoms_input = AtomsInput(atoms=self.sc_atoms)

    def create_max_val(self):
        dryrun = DryrunJob(atoms=self.atoms_input, scheduler=self.scheduler_dryrun, job_desc=self.job_single_node, is_spinorbit=self.scf_is_spinorbit)
        dryrun.create()
        dryrun.run(0.0)
        self.max_val = dryrun.get_max_val()
        dryrun.remove()

    def create_kpath(self):
        self.kpath_obj = KPath(
            atoms=self.uc_atoms,
            path_special_points=self.path_special_points,
            path_segment_npoints=self.path_segment_npoints,
            path_total_npoints=self.path_total_npoints,
        )
        save_obj(self.kpath_obj, 'bandpath.pkl')
        # self.Kpath, self.Gpath = self.kpath_obj.get_sc_path(self.sc_grid)

    def create_jobs_input(self, save=True):
        self.relax = RelaxInput(
            max_val=self.max_val,
            job_desc=self.job_para,
            relax_type=self.relax_type,
            extra_control_args=self.relax_extra_control_args,
            extra_system_args=self.relax_extra_system_args,
            extra_electron_args=self.relax_extra_electrons_args,
        )

        self.scf = ScfInput(
            kdim=self.scf_kgrid,
            ecutwfc=self.scf_cutoff,
            job_desc=self.job_para,
            is_spinorbit=self.scf_is_spinorbit,
            num_val_bands=self.max_val,
            extra_control_args=self.scf_extra_control_args,
            extra_system_args=self.scf_extra_system_args,
            extra_electron_args=self.scf_extra_electrons_args,
        )

        self.dfpt = DfptInput(
            atoms=self.atoms_input,
            qgrid=self.dfpt_qgrid,
            conv_threshold=self.dfpt_conv_threshold,
            job_desc=self.job_big_para_k,
            extra_args=self.dfpt_extra_args,
        )

        self.phbands = PhbandsInput(
            kpath=self.kpath_obj,
            job_desc=self.job_single_node,
            extra_q2r_args=self.phbands_extra_q2r_args,
            extra_matdyn_args=self.phbands_extra_matdyn_args,
        )

        self.phdos = PhdosInput(
            qdim=self.dos_kdim,
            job_desc=self.job_single_node,
            extra_q2r_args=self.phdos_extra_q2r_args,
            extra_matdyn_args=self.phdos_extra_matdyn_args,
        )

        self.phmodes = PhmodesInput(
            qidx=self.phmodes_qpt_idx,
            job_desc=self.job_single_node,
            extra_args=self.phmodes_extra_args,
        )

        self.dos = DosInput(
            kdim=self.dos_kdim,
            bands=self.dftelbands_cond + self.max_val,
            job_desc=self.job_big_para,
            extra_control_args=self.dos_extra_control_args,
            extra_system_args=self.dos_extra_system_args,
            extra_electrons_args=self.dos_extra_electrons_args,
            extra_dos_args=self.dos_extra_args,
            extra_pdos_args=self.pdos_extra_args,
        )

        self.dftelbands = DftelbandsInput(
            kpath=self.kpath_obj,
            nbands=self.dftelbands_cond + self.max_val,
            job_desc=self.job_para,
            job_pw2bgw_desc=self.job_single_node,
            extra_control_args=self.dftelbands_extra_control_args,
            extra_system_args=self.dftelbands_extra_system_args,
            extra_electron_args=self.dftelbands_extra_electrons_args,
        )

        self.kpdos = KpdosInput(
            job_desc = self.job_big_para,
            extra_args=self.kpdos_extra_args,
        )

        self.wannier = WannierInput(
            atoms=self.atoms_input,
            kdim=self.wannier_kdim,
            num_bands=self.wannier_bands_cond + self.max_val,
            num_wann=self.wannier_bands_cond + self.max_val,
            job_wfnwan_desc=self.job_para,
            job_pw2wan_desc=self.job_single_node,
            job_wan_desc=self.job_para_epwk,
            extra_control_args=self.wannier_extra_control_args,
            extra_system_args=self.wannier_extra_system_args,
            extra_electron_args=self.wannier_extra_electrons_args,
            extra_args=self.wannier_extra_args,
        )

        self.wfn = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfn_qe_kdim,
            qshift=(0.0, 0.0, 0.0),
            is_reduced=self.wfn_qe_sym,
            bands=self.wfn_qe_cond + self.max_val,
            job_wfn_desc=self.job_para_k,
            job_pw2bgw_desc=self.job_single_node,
            job_parabands_desc=self.job_big_para,
            parabands_bands=self.wfn_para_cond + self.max_val,
            extra_control_args=self.wfn_extra_control_args,
            extra_system_args=self.wfn_extra_system_args,
            extra_electrons_args=self.wfn_extra_electrons_args,
            extra_parabands_args=self.wfn_extra_parabands_args,
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
            job_desc=self.job_para_epwk,
            skipped_bands=skipped_bands,     # The input bands are 1 to 14, which are fine.
            extra_args=self.epw_extra_args,
        )

        self.wfnq = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnq_qe_kdim,
            qshift=self.qshift,
            is_reduced=self.wfnq_qe_sym,
            bands=self.wfnq_qe_cond + self.max_val,
            job_wfn_desc=self.job_para_k,
            job_pw2bgw_desc=self.job_single_node,
            extra_control_args=self.wfnq_extra_control_args,
            extra_system_args=self.wfnq_extra_system_args,
            extra_electrons_args=self.wfnq_extra_electrons_args,
        )

        self.wfnfi = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnfi_qe_kdim,
            qshift=(0.0, 0.0, 0.000),
            is_reduced=self.wfnfi_qe_sym,
            bands=self.wfnfi_qe_cond,
            job_wfn_desc=self.job_para_k,
            job_pw2bgw_desc=self.job_single_node,
            extra_control_args=self.wfnfi_extra_control_args,
            extra_system_args=self.wfnfi_extra_system_args,
            extra_electrons_args=self.wfnfi_extra_electrons_args,
        )

        self.wfnqfi = WfnGeneralInput(
            atoms=self.atoms_input,
            kdim=self.wfnqfi_qe_kdim,
            qshift=(0.0, 0.0, 0.001),
            is_reduced=self.wfnqfi_qe_sym,
            bands=self.wfnqfi_qe_cond,
            job_wfn_desc=self.job_para_k,
            job_pw2bgw_desc=self.job_single_node,
            extra_control_args=self.wfnqfi_extra_control_args,
            extra_system_args=self.wfnqfi_extra_system_args,
            extra_electrons_args=self.wfnqfi_extra_electrons_args,
        )

        self.epsilon = EpsilonInput(
            bands=self.epssig_bands_cond + self.max_val,
            cutoff=self.epssig_cutoff,
            wfn_link=self.epssig_wfnlink,
            wfnq_link=self.epssig_wfnqlink,
            job_desc=self.job_para,
            extra_args=self.eps_extra_args,
        )

        self.sigma = SigmaInput(
            bands=self.epssig_bands_cond + self.max_val,
            band_min=self.max_val - self.sig_band_val + 1,
            band_max=self.max_val + self.sig_band_cond,
            cutoff=self.epssig_cutoff,
            wfn_inner_link=self.epssig_wfnlink,
            job_desc=self.job_para,
            extra_args=self.sig_extra_args,
        )

        self.inteqp = InteqpInput(
            val_bands_coarse=self.inteqp_band_val,
            cond_bands_coarse=self.dftelbands_cond-1,
            val_bands_fine=self.inteqp_band_val,
            cond_bands_fine=self.dftelbands_cond-1,
            wfn_co_link=self.inteqp_wfn_co_link,
            wfn_fi_link=self.inteqp_wfn_fi_link,
            job_desc=self.job_para,
            extra_args=self.inteqp_extra_args,
        )

        self.kernel = KernelInput(
            val_bands_coarse=self.abs_val_bands,
            cond_bands_coarse=self.abs_cond_bands,
            Qshift=(0.0, 0.0, 0.0),
            wfn_co_link=self.abs_wfn_co_link,
            wfnq_co_link=self.abs_wfnq_co_link,
            job_desc=self.job_para,
            extra_args=self.ker_extra_args,
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
            pol_dir=self.abs_pol_dir,
            job_desc=self.job_para,
            extra_args=self.abs_extra_args,
        )
        
        self.plotxct = PlotxctInput(
            hole_position=self.plotxct_hole,
            supercell_size=self.plotxct_sc,
            state=self.plotxct_state,
            wfn_fi_link=self.plotxct_wfn_fi_link,
            wfnq_fi_link=self.plotxct_wfnq_fi_link,
            job_desc=self.job_para,
            extra_args=self.plotxct_extra_args,
        )

        self.bseq = BseqInput(
            atoms=self.atoms_input,
            Qdim=self.bseq_Qdim,
            job_desc=self.job_big_para,
        )

        self.xctphbgw = XctPhInput(
            job_desc=self.job_single_task,
            epw_qgrid=self.bseq_Qdim,
            num_epw_val_bands=self.abs_val_bands,
            num_epw_cond_bands=self.abs_cond_bands,
            num_exciton_states=self.abs_nevec,
        )

        self.esf = EsfInput(
            job_desc=self.job_single_task,
        )

        self.xctpol = XctPolInput(
            max_error=self.xctpol_max_error,
            max_steps=self.xctpol_max_steps,
            job_desc=self.job_single_task,
        )

        self.input: Input = Input(
            scheduler=self.scheduler,
            atoms=self.atoms_input,
            relax=self.relax,
            scf=self.scf,
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
            xctph=self.xctphbgw,
            esf=self.esf,
            xctpol=self.xctpol,
        )
        if save: save_obj(self.input, 'input.pkl')

    def create_input(self, save=True):
        
        self.create_atoms()

        if not self.atoms_skip_pseudo_generation: self.create_pseudos()

        self.create_kpath()

        self.create_atoms_input()

        self.create_max_val()

        self.create_jobs_input(save)

    def get_flowmanage(self, list_of_step_classes: list, save_pkl: bool =True) -> FlowManage:
        self.create_input(save_pkl)

        list_of_steps = [step_class(self.input) for step_class in list_of_step_classes]
        self.flowmanage: FlowManage = FlowManage(list_of_steps)
        if save_pkl: save_obj(self.flowmanage, 'flowmanage.pkl'); save_obj(self, 'fullgridflow.pkl')
        return self.flowmanage

#endregion
