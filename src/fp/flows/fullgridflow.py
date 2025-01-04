#region: Modules.
from ase import Atoms 
from ase.io import write, read
import numpy as np 
from ase.build import make_supercell
import yaml 
from typing import List
from fp.inputs.atoms import AtomsInput
from fp.inputs.relax import RelaxInput
from fp.inputs.scf import ScfInput
from fp.inputs.phbands import PhbandsInput
from fp.inputs.dftelbands import DftelbandsInput
from fp.inputs.wannier import WannierInput
from fp.inputs.wfngeneral import WfnGeneralInput
from fp.inputs.epw import EpwInput
from fp.inputs.epsilon import EpsilonInput
from fp.inputs.sigma import SigmaInput
from fp.inputs.abs import PlotxctInput
from fp.inputs.bseq import BseqInput
from fp.inputs.input_main import Input
from fp.jobs.dryrun import DryrunJob
from fp.flows.flow_manage import FlowManage
from fp.io.pkl import save_obj, load_obj
from fp.structure.kpath import KPath
from fp.schedulers.scheduler import Scheduler
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
        self.input_dict: dict = None
        self.input: Input = None
    
        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def from_yml(filename):
        '''
        Generate a fullgrid flow object from a yml file.
        '''
        # Open and read the YAML file
        with open(filename, 'r') as file:
            yml_data: dict = yaml.safe_load(file)

        fullgridflow: FullGridFlow = FullGridFlow(input_dict=yml_data)

        return fullgridflow

    def create_interactive(self):
        Scheduler.from_input_dict(self.input_dict).create_interactive()

    def create_pseudos(self):
        FlowManage.create_pseudos(
            atoms=self.input.atoms.atoms, 
            is_fr=self.input_dict['scf']['is_spinorbit'], 
            xc_type=self.input_dict['scf']['xc_type'], 
            pseudos_dict=self.input_dict['atoms']['pseudos'],
        )

    def create_max_val(self):
        dryrun = DryrunJob(input=self.input)
        dryrun.create()
        dryrun.run(0.0)
        self.input_dict['total_valence_bands'] = dryrun.get_max_val()
        dryrun.remove()

    def create_kpath(self):
        self.kpath_obj = KPath(self.input_dict, self.input.atoms.atoms)
        save_obj(self.kpath_obj, 'bandpath.pkl')
        # self.Kpath, self.Gpath = self.kpath_obj.get_sc_path(self.sc_grid)

    def create_input(self):
        # Components. 
        atoms = AtomsInput(self.input_dict)
        write('flow_struct.xsf', atoms.atoms)
        relax = RelaxInput(self.input_dict)
        scf = ScfInput(self.input_dict)
        phbands = PhbandsInput(
            input_dict=self.input_dict, 
            atoms=atoms.atoms
        )
        dftelbands = DftelbandsInput(self.input_dict)
        wannier = WannierInput(
            input_dict=self.input_dict,
            atoms_input=atoms,
        )
        wfn = WfnGeneralInput(
            input_dict=self.input_dict,
            atoms_input=atoms,
            wfn_type='wfn',
        )
        epw = EpwInput(self.input_dict)
        wfnq = WfnGeneralInput(
            input_dict=self.input_dict,
            atoms_input=atoms,
            wfn_type='wfnq',
        )
        wfnfi = WfnGeneralInput(
            input_dict=self.input_dict,
            atoms_input=atoms,
            wfn_type='wfnfi',
        )
        wfnqfi = WfnGeneralInput(
            input_dict=self.input_dict,
            atoms_input=atoms,
            wfn_type='wfnqfi',
        )
        eps = EpsilonInput(self.input_dict)
        sig = SigmaInput(self.input_dict)
        plotxct = PlotxctInput(self.input_dict)
        bseq = BseqInput(self.input_dict)

        # Aggregate.
        self.input: Input = Input(
            input_dict=self.input_dict,
            atoms=atoms,
            relax=relax,
            scf=scf,
            phbands=phbands,
            dftelbands=dftelbands,
            wannier=wannier,
            wfn=wfn,
            epw=epw,
            wfnq=wfnq,
            wfnfi=wfnfi,
            wfnqfi=wfnqfi,
            eps=eps,
            sig=sig,
            plotxct=plotxct,
            bseq=bseq,
        )
        save_obj(self.input_dict, 'input_dict.pkl')
        save_obj(self.input, 'input.pkl')

    def build_steps(self):
        self.create_input()
        self.create_pseudos()
        self.create_kpath()
        self.create_max_val()

    def get_flowmanage(self, list_of_step_classes: list) -> FlowManage:
        self.build_steps()

        list_of_steps = [step_class(self.input) for step_class in list_of_step_classes]
        self.flowmanage: FlowManage = FlowManage(list_of_steps)
        save_obj(self.flowmanage, 'flowmanage.pkl')
        save_obj(self, 'fullgridflow.pkl')
        
        return self.flowmanage

#endregion
