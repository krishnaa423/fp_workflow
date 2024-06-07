#region: Modules.
#endregion

#region: Variables.
#endregion

#region: Functions.
def update():
    pass

def run():
    pass 

def save():
    pass 

def step_single(parameter):
    pass 

def main():
    conv_params = {
        'val_coarse': [],
        'cond_coarse': [],
        'val_fine': [],
        'cond_fine': [],
    }
    
    for key in conv_params.keys():
        for value in conv_params[key]:
            step_single({key: value})
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion