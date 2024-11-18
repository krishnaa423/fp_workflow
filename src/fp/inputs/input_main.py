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
from fp.inputs.inteqp import *
from fp.inputs.kernel import *
from fp.inputs.abs import *
from fp.inputs.bseq import *
from fp.inputs.esf import *
from fp.inputs.xctph import *
from fp.inputs.xctpol import *

from fp.schedulers import *
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

        relax: RelaxInput,
        scf: ScfInput,

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
        inteqp: InteqpInput,
        kernel: KernelInput,
        absorption: AbsorptionInput,
        plotxct: PlotxctInput,
        bseq: BseqInput,
        
        xctph: XctPhInput,
        esf: EsfInput,
        xctpol: XctPolInput,
    ):
        self.scheduler: Scheduler = scheduler
        self.atoms: AtomsInput = atoms 

        self.relax: RelaxInput = relax 
        self.scf: ScfInput = scf

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
        self.inteqp: InteqpInput = inteqp
        self.kernel: KernelInput = kernel
        self.absorption: AbsorptionInput = absorption
        self.plotxct: PlotxctInput = plotxct
        self.bseq: BseqInput = bseq
        
        self.xctph: XctPhInput = xctph
        self.esf: EsfInput = esf

        # self.pol: PolInput = pol
        self.xctpol: XctPolInput = xctpol
#endregion