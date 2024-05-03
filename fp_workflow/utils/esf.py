#region: Modules.
import h5py 
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    pass  
#endregion

#region: Classes.
class DftStruct:
    def get_struct(self):
        pass 

class DftForce:
    def get_dftforce(self):
        pass 

class Eqp:
    def get_dfteig(self):
        pass 
    
    def get_gweig(self):
        pass 

class Elph:
    def get_elph(self):
        pass

class Xct:
    def get_xcteig(self):
        pass 
    
    def get_xctevec(self):
        pass 
    
class Esf:
    def assemble_components(self):
        pass 
    
    def do_calculation(self):
        pass 
    
    def write(self):
        
        with h5py.File('esf.h5', 'w') as f:
            ds_positions = f.create_dataset('positions')
            ds_dft_force = f.create_dataset('dft_force')
            ds_excited_force = f.create_dataset('excited_force')
            
            # ds_positions = 
            # ds_dft_force = 
            # ds_excited_force = 
    
    def get_esf(self):
        self.assemble_components()
        self.do_calculation()
        self.write()
#endregion

#region: Main.
if __name__ == '__main__':
    main()
#endregion
