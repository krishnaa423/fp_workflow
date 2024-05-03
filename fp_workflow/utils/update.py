#region: Modules.
from ase import Atoms 
import workflow
#endregion

#region: Variables.
#endregion

#region: Functions.
def read_step():
    pass 

def get_struct_from_step():
    step = read_step()
    pass 


def main():
    
    # Create struct. 
    struct = get_struct_from_step()
    
    # Remove current and create new workflow. 
    workflow.remove_workflow()
    workflow.create_workflow(struct=struct)      # pass in struct as a parameter. Both struct and input json will be optionals with defaults. 
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__ == "__main__":
    main()
#endregion
