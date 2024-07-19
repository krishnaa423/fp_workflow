#region: Modules.
from fp.flows import *
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    flow: FlowManage = FlowManage.load_flow('./flowmanage.pkl')
    flow.run()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion