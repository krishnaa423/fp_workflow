#region: Modules.
from fp.inputs.atoms import *
from fp.inputs.scf import *
from fp.inputs.relax import *
from fp.inputs.dfpt import *
from fp.inputs.phbands import *
from fp.inputs.phdos import *
from fp.inputs.phmodes import *
from fp.inputs.dos import *
from fp.inputs.dftelbands import *
from fp.inputs.kpdos import *
from fp.inputs.wannier import *
from fp.inputs.wfngeneral import *
from fp.inputs.epw import *
from fp.inputs.epsilon import *
from fp.inputs.sigma import *

from fp.schedulers import *
import pickle
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Input:
    def __init__(
        self,
        scheduler: Scheduler,
        atoms: AtomsInput,
        scf: ScfInput,
        relax: RelaxInput,
        dfpt: DfptInput,
        phbands: PhbandsInput,
        phdos: PhdosInput,
        phmodes: PhmodesInput,
        dos: DosInput,
        dftelbands: DftelbandsInput,
        kpdos: KpdosInput,
        wannier: WannierInput,
        wfn: WfnGeneralInput,
        epw: EpwInput,
        wfnq: WfnGeneralInput,
        wfnfi: WfnGeneralInput,
        wfnqfi: WfnGeneralInput,
        epsilon: EpsilonInput,
        sigma: SigmaInput,
    ):
        self.scheduler: Scheduler = scheduler
        self.atoms: AtomsInput = atoms 
        self.scf: ScfInput = scf
        self.relax: RelaxInput = relax 
        self.dfpt: DfptInput = dfpt
        self.phbands: PhbandsInput = phbands
        self.phdos: PhdosInput = phdos
        self.phmodes: PhmodesInput = phmodes

        self.dos: DosInput = dos
        self.dftelbands: DftelbandsInput = dftelbands
        self.kpdos: KpdosInput = kpdos
        self.wannier: WannierInput = wannier
        
        self.wfn: WfnGeneralInput = wfn
        self.epw: EpwInput = epw
        self.wfnq: WfnGeneralInput = wfnq
        self.wfnfi: WfnGeneralInput = wfnfi
        self.wfnqfi: WfnGeneralInput = wfnqfi

        self.epsilon: EpsilonInput = epsilon
        self.sigma: SigmaInput = sigma
        # self.kernel: KernelInput = kernel
        # self.absorption: AbsorptionInput = absorption
        # self.plotxct: PlotxctInput = plotxct
        # self.bseq: BseQInput = bseq
        # self.xctph: XctphInput = xctph

        # self.esf: EsfInput = esf
        # self.esd: EsdInput = esd

        # self.pol: PolInput = pol
        # self.xctpol: XctpolInput = xctpol

    def save(self, filename='input.pkl'):
        with open(filename, 'wb') as f: pickle.dump(self, f)
#endregion