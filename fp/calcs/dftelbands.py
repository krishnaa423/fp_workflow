#region: Modules.
from fp.inputs import *
from fp.io import *
from fp.flows import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Dftelbands:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_dftelbands = \
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
nbnd={self.input.dftelbands.nbands}
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

{self.input.dftelbands.get_kgrid_str()}
'''
        
        self.job_dftelbands = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dftelbands.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dftelbands.job_desc)}pw.x < dftelbands.in &> dftelbands.in.out  
cp ./tmp/struct.xml ./dftelbands.xml 
'''
        
        self.input_dftelbands_pw2bgw = \
f'''&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFN_dftelbands'
wfng_kgrid=.true.
wfng_nk1=0
wfng_nk2=0
wfng_nk3=0
wfng_dk1=0.0
wfng_dk2=0.0
wfng_dk3=0.0
/
'''
        
        self.job_dftelbands_pw2bgw = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.dftelbands.job_pw2bgw_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.dftelbands.job_pw2bgw_desc)}pw2bgw.x -pd .true. < dftelbands_pw2bgw.in &> dftelbands_pw2bgw.in.out 
cp ./tmp/WFN_dftelbands ./
wfn2hdf.x BIN WFN_dftelbands WFN_dftelbands.h5
'''
        
        self.jobs = [
            'job_dftelbands.sh',
            'job_dftelbands_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f('dftelbands.in', self.input_dftelbands)
        write_str_2_f('job_dftelbands.sh', self.job_dftelbands)
        write_str_2_f('dftelbands_pw2bgw.in', self.input_dftelbands_pw2bgw)
        write_str_2_f('job_dftelbands_pw2bgw.sh', self.job_dftelbands_pw2bgw)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_dftelbands.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_dftelbands_pw2bgw.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'dftelbands.in*',
            'job_dftelbands.sh',
            'dftelbands_pw2bgw.in*',
            'job_dftelbands_pw2bgw.sh',
            'WFN_dftelbands',
            'WFN_dftelbands.h5',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf dftelbands.in')
        os.system('rm -rf job_dftelbands.sh')
        os.system('rm -rf dftelbands_pw2bgw.in')
        os.system('rm -rf job_dftelbands_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf dftelbands.in.out')
        os.system('rm -rf dftelbands_pw2bgw.in.out')
        os.system('rm -rf dftelbands.xml')
        os.system('rm -rf kgrid.inp kgrid.log kgrid.out')
        os.system('rm -rf WFN_dftelbands')
        os.system('rm -rf WFN_dftelbands.h5')
        
        os.system('rm -rf uc_kpath.txt')
        os.system('rm -rf sc_grid.txt')
        os.system('rm -rf sc_Kpath.txt')
        os.system('rm -rf sc_Gshift.txt')
#endregion
