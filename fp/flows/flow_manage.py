#region: Modules.
from ase import Atoms 
from pkg_resources import resource_filename
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class FlowManage:
    def __init__(
        self,
        list_of_steps,
        is_fr=False,    # Is fully relativistic for the pseudopotentials?
    ):
        self.list_of_steps: list = list_of_steps
        self.is_fr: bool = is_fr

    @staticmethod
    def create_pseudos(atoms: Atoms, is_fr: bool=False):
        sr_fr_str = 'fr' if is_fr else 'sr'

        pkg_dir = resource_filename('fp', '')
        pseudo_dir = pkg_dir + '/pseudos/ONCVPSP/sg15'

        os.system('mkdir -p ./ONCVPSP')
        os.system('mkdir -p ./ONCVPSP/sg15')

        symbols = atoms.get_chemical_symbols()

        for sym in symbols:
            source_file = pseudo_dir + f'/{sym}_ONCV_PBE_{sr_fr_str}.upf'
            dest_file = './ONCVPSP/sg15' + f'/{sym}_ONCV_PBE_{sr_fr_str}.upf'
            os.system(f'cp {source_file} {dest_file}')

    def create(self):
        assert len(self.list_of_steps)>=1, 'Number of steps/jobs should be greater than 0.'

        # Add the pseudopotential files.
        self.create_pseudos(self.list_of_steps[0].input.atoms.atoms, self.is_fr)

        for step in self.list_of_steps:
            step.create()

    def run(self, total_time=0.0):
        total_time: float = total_time

        for step in self.list_of_steps:
            total_time = step.run(total_time)

        # Write the total workflow run time. 
        print(f'Done whole worflow in {total_time:15.10f} seconds.\n\n', flush=True)

    def save(self, folder):       
        for step in self.list_of_steps:
            step.save(folder)

    def remove(self):
        for step in self.list_of_steps:
            step.remove()
         
#endregion