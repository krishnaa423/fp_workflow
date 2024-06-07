#region: Modules.
import subprocess
import os 
import time
from ..workflow import Scf, DftElbands
import pickle 
import argparse 
#endregion

#region: Variables.
#endregion

#region: Functions.
def make_dir(name):
    if not os.path.exists(name): os.system(f'mkdir {name}')

def main():
    
    parser = argparse.ArgumentParser(description='DFT convergence.')
    parser.add_argument('--run', action='store_true', help='Run scripts.')
    parser.add_argument('--plot', action='store_true', help='Generate plot. ')
    
    
    args = parser.parse_args()
    
    ecuts = [10.0, 15.0, 20.0]
    dft_conv = DftConvergence(
        ecuts=ecuts
    )
        
    if args.run:
        dft_conv.run_all()
        
    if args.plot:
        dft_conv.gen_plot()
#endregion

#region: Classes.
class DftConvergence:
    def __init__(
        self,
        ecuts,
    ):
        self.ecuts = ecuts
        
        with open('../input.pkl', 'rb') as rb: self.input = pickle.load(rb)
        self.list_of_steps = [
            Scf(input=self.input),
            DftElbands(input=self.input),
        ]
    
    def update_single(self, ecut):
        pattern = r'ecutwfc=.(\<>*\n)'
    
    def run_single(self):
        os.system('cd ..')
        for step in self.list_of_steps:
            step.run()
        os.system('cd ./utils')
    
    def save_single(self, ecut):
        os.system('cd ..')
        folder_name = f'./convergence/dft/ecut{ecut:10.5f}'
        make_dir('./convergence')
        make_dir('./convergence/dft')
        make_dir(folder_name)
        os.system('cd ./utils')
        
        # Save all the necessary folders.
        for step in self.list_of_steps:
            step.save(folder_name)
        
    def run_all(self, ecut):
        for ecut in self.ecuts:
            self.update_single(ecut=ecut)
            self.run_single()
            self.save_single(ecut=ecut)
    
    def get_energies(self):
        pass 
    
    def gen_plot(self):
        pass 
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion