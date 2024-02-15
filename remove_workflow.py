import os

inodes_not_to_remove = [
    'abinit_pseudos',
    '.vscode',
    'qe_pseudos',
    'analysis',
    'backup',
    'create_workflow.py',
    'remove_workflow.py',
    'run_workflow.py',
    'abi2bgw_nc_h5.py',
    'gen_struct.py',
    'esf.py',
    'step.py',
    'update.py',
    'input.json',
    'gen_input.json',
    'step.h5',
    'job_run_workflow.sh',
    'job_interactive.sh',
    'sample.py',
    'explore.py',
    'test.py',
    'test.abi',
    'mrgddb.files',
    'anaddb.files',
    'anaddb.abi',
    'wannier90.win',
    'test_gpaw.py',
    'test_nn.py',
    'test_api.py',
]

inodes = os.listdir('./')

for inode in inodes:
    if not inode in inodes_not_to_remove:
        os.system(f'rm -rf {inode}')