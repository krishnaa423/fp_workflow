#region: Modules.
import h5py 
import numpy as np 
#endregion

#region: Variables.
#endregion

#region: Functions.
def get_iter():
    pass 

def main():
    # Get inputs. 
    esf = Esf()
    iter = get_iter()
    
    # Bfgs. 
    bfgs = Bfgs(esf, iter)
    bfgs.assemble()
    bfgs.run_step()
    bfgs.write_step()
#endregion

#region: Classes.
class Bfgs:
    def __init__(self, esf, iter):
        self.esf = esf
        self.iter = iter 
        self.x: np.ndarray = np.zeros()
        self.x_next: np.ndarray = np.zeros()
        self.f_prev: np.ndarray = np.zeros()
        self.f: np.ndarray = np.zeros()
        self.H_next: np.ndarray = np.zeros()
        self.H: np.ndarray = np.zeros()
        self.p: np.ndarray = np.zeros()
        self.alpha: np.ndarray = np.zeros()
        self.s: np.ndarray = np.zeros()
        self.y: np.ndarray = np.zeros()
        
        # Solver params. 
        self.max_iter = 0
        self.max_error = 1.0e-1
        
        self.read_input()
    
    def read_input(self):
        pass 
    
    def assemble(self):
        pass 
    
    def write_step(self):
        with h5py.File('step.h5', 'w') as f:
            ds_x = f.create_dataset('x')
            ds_f = f.create_dataset('f')
            ds_H = f.create_dataset('H')
            ds_p = f.create_dataset('p')
            ds_alpha = f.create_dataset('alpha')
            ds_s = f.create_dataset('s')
            ds_y = f.create_dataset('y')
            ds_y = f.create_dataset('iter')
            
            # ds_x = 
            # ds_f = 
            # ds_H = 
            # ds_p = 
            # ds_alpha = 
            # ds_s = 
            # ds_y = 
            # ds_iter = 
        
    def run_step(self):
        self.p = - np.matmul(self.H, self.x)
        self.s = self.alpha*self.p
        self.x_next = self.x + self.s 
        self.y = self.f - self.f_prev
        self.H_next = self.H  # + two terms. 
        
        self.iter += 1
        
class Esf:
    def get_positions(self):
        pass 
    
    def get_force(self):
        pass 
#endregion

#region: Main.
if __name__ == '__main__':
    pass 
#endregion