#region: Modules.
from fp.flows import *
#endregion

#region: Variables.
#endregion

#region: Functions.
def main():
    flow: FlowManage = FlowManage.load_flow('./flowmanage.pkl')
    flow.remove(pkl_too=False)
    os.system('rm -rf job_all*')
    os.system('rm -rf ./test_save_folder')
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion