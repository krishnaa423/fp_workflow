#region: Modules.
from fp.flows import *
from fp.io import *
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    flow: FlowManage = load_obj('flowmanage.pkl')
    flow.run()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion