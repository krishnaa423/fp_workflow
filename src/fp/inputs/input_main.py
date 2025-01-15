#region modules
from fp.inputs.abs import PlotxctInput
from fp.inputs.atoms import AtomsInput
from fp.inputs.bseq import BseqInput
from fp.inputs.dftelbands import DftelbandsInput
from fp.inputs.epsilon import EpsilonInput
from fp.inputs.epw import EpwInput
from fp.inputs.phbands import PhbandsInput
from fp.inputs.relax import RelaxInput
from fp.inputs.scf import ScfInput
from fp.inputs.sigma import SigmaInput
from fp.inputs.wfngeneral import WfnGeneralInput
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class Input:
    def __init__(
        self,
        input_dict: dict=None,
        atoms: AtomsInput=None,
        relax: RelaxInput=None,
        scf: ScfInput=None,
        phbands: PhbandsInput=None,
        dftelbands: DftelbandsInput=None,
        wfn: WfnGeneralInput=None,
        epw: EpwInput=None,
        wfnq: WfnGeneralInput=None,
        wfnfi: WfnGeneralInput=None,
        wfnqfi: WfnGeneralInput=None,
        eps: EpsilonInput=None,
        sig: SigmaInput=None,
        plotxct: PlotxctInput=None,
        bseq: BseqInput=None,
        **kwargs,
    ):
        self.input_dict: dict = input_dict
        self.atoms: AtomsInput = atoms
        self.relax: RelaxInput = relax
        self.scf: ScfInput = scf
        self.phbands: PhbandsInput = phbands
        self.dftelbands: DftelbandsInput = dftelbands
        self.wfn: WfnGeneralInput = wfn
        self.epw: EpwInput = epw
        self.wfnq: WfnGeneralInput = wfnq
        self.wfnfi: WfnGeneralInput = wfnfi
        self.wfnqfi: WfnGeneralInput = wfnqfi
        self.eps: EpsilonInput = eps
        self.sig: SigmaInput = sig
        self.plotxct: PlotxctInput = plotxct
        self.bseq: BseqInput = bseq

        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def from_dict(input_dict: dict):
        # Components. 
        atoms = AtomsInput(input_dict)
        
        relax = RelaxInput(input_dict)
        scf = ScfInput(input_dict)
        phbands = PhbandsInput(
            input_dict=input_dict, 
            atoms=atoms.atoms
        )
        dftelbands = DftelbandsInput(input_dict)
        wfn = WfnGeneralInput(
            input_dict=input_dict,
            atoms_input=atoms,
            wfn_type='wfn',
        )
        epw = EpwInput(input_dict)
        wfnq = WfnGeneralInput(
            input_dict=input_dict,
            atoms_input=atoms,
            wfn_type='wfnq',
        )
        wfnfi = WfnGeneralInput(
            input_dict=input_dict,
            atoms_input=atoms,
            wfn_type='wfnfi',
        )
        wfnqfi = WfnGeneralInput(
            input_dict=input_dict,
            atoms_input=atoms,
            wfn_type='wfnqfi',
        )
        eps = EpsilonInput(input_dict)
        sig = SigmaInput(input_dict)
        plotxct = PlotxctInput(input_dict)
        bseq = BseqInput(input_dict)

        input = Input(
            input_dict=input_dict,
            atoms=atoms,
            relax=relax,
            scf=scf,
            phbands=phbands,
            dftelbands=dftelbands,
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

        return input
    
    def update_qe_args_dict(self, args_dict: dict, args_type: str, qedict_to_update: dict):
        if args_dict is not None:
            if args_type=='override':
                qedict_to_update = args_dict.copy()
            if args_type=='extra':
                if 'namelists' in args_dict:
                    for key, value in args_dict['namelists'].items():
                        if qedict_to_update.get('namelists') is None: qedict_to_update['namelists'] = {}
                        if qedict_to_update.get('namelists', {}).get(key) is None: qedict_to_update['namelists'][key] = {}
                        qedict_to_update['namelists'][key].update(args_dict['namelists'][key])
                if 'blocks' in args_dict:
                    if qedict_to_update.get('blocks') is None: qedict_to_update['blocks'] = {}
                    qedict_to_update['blocks'].update(args_dict['blocks'])

        return qedict_to_update

#endregions
