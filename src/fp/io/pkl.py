#region: Modules.
import dill as pickle
#endregion

#region: Variables.
#endregion

#region: Functions.
def save_obj(obj, filename):

    with open(filename, 'wb') as w: 
        pickle.dump(obj, w)

def load_obj(filename):
    with open(filename, 'rb') as w: 
        obj = pickle.load(w)

    return obj

#endregion

#region: Classes.
#endregion
