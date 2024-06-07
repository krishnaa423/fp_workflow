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
        'screened_cutoff': [],
        'sigma_bands': [],
        'epsilon_bands': [],
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