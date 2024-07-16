#region: Modules.
from fp.inputs.atoms import *
from fp.inputs.scf import *
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
    ):
        self.scheduler: Scheduler = scheduler
        self.atoms: AtomsInput = atoms 
        self.scf: ScfInput = scf
#endregion