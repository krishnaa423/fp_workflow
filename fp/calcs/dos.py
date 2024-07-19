#region: Modules.
from fp.inputs import *
from fp.io import *
from fp.flows import *
from pkg_resources import resource_filename
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Dos:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_wfndos = \
f'''&CONTROL
outdir='./tmp'
prefix='struct'
pseudo_dir='./ONCVPSP/sg15'
calculation='bands'
tprnfor=.true. 
/

&SYSTEM
ibrav=0
ntyp={self.input.atoms.get_ntyp()}
nat={self.input.atoms.get_nat()}
nbnd={self.input.dos.bands}
ecutwfc={self.input.scf.ecutwfc}
!noncolin=.true.
!lspinorb=.true. 
/

&ELECTRONS
/

&CELL
/

&IONS
/

CELL_PARAMETERS angstrom
{self.input.atoms.get_scf_cell()}

ATOMIC_SPECIES
{self.input.atoms.get_scf_atomic_species()}

ATOMIC_POSITIONS angstrom 
{self.input.atoms.get_scf_atomic_positions()}

K_POINTS automatic 
{int(self.input.dos.kdim[0])} {int(self.input.dos.kdim[1])} {int(self.input.dos.kdim[2])} 0 0 0
'''
        
        self.job_wfndos = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dos.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dos.job_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.dos.job_desc)} < wfndos.in &> wfndos.in.out 
'''
        
        self.input_dos = \
f'''&DOS
outdir='./tmp'
prefix='struct'
fildos='struct_dos.dat'
/
'''
        
        self.job_dos = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dos.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dos.job_desc)}dos.x -pd .true. < dos.in &> dos.in.out 
'''
        self.jobs = [
            'job_wfndos.sh',
            'job_dos.sh'
        ]

    def create(self):
        write_str_2_f('wfndos.in', self.input_wfndos)
        write_str_2_f('job_wfndos.sh', self.job_wfndos)
        write_str_2_f('dos.in', self.input_dos)
        write_str_2_f('job_dos.sh', self.job_dos)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_wfndos.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_dos.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wfndos.in*',
            'dos.in*',
            'struct_dos.dat',
            'job_dos.sh',
            'job_wfndos.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wfndos.in')
        os.system('rm -rf dos.in')
        os.system('rm -rf job_dos.sh')
        os.system('rm -rf job_wfndos.sh')
        
        os.system('rm -rf struct_dos.dat')
        os.system('rm -rf dos.in.out')
        os.system('rm -rf wfndos.in.out')
#endregion
