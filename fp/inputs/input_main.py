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
        scf: ScfInput,
        relax: RelaxInput,
        dfpt: DfptInput,
        phbands: PhbandsInput,
        phdos: PhdosInput,
        phmodes: PhmodesInput,
        dos: DosInput,
        dftelbands: DftelbandsInput,
        kpdos: KpdosInput,
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
#endregion