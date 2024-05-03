#region: Modules.
import workflow 
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
def read_iter_from_step():
    pass 

def main():
    # Get folder name. 
    iter = read_iter_from_step()
    folder_name = f'./savestep/iter{iter}'
    
    # Remove current and create new folder. 
    os.system(f'rm -rf {folder_name}')
    os.system(f'mkdir {folder_name}')
    
    # Save stuff to the folder. 
    workflow.save_step(folder_name)     
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__ == "__main__":
    main()
#endregion