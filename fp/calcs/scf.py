#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import *
from fp.flows.run import *
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Scf:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_scf: str = \
f'''&CONTROL
outdir='./tmp'
prefix='struct'
pseudo_dir='./ONCVPSP/sg15'
calculation='scf'
tprnfor=.true.
/

&SYSTEM
ibrav=0
ntyp={self.input.atoms.get_ntyp()}
nat={self.input.atoms.get_nat()}
ecutwfc={self.input.scf.ecutwfc}
!noncolin=.true.
!lspinorb=.true. 
/

&ELECTRONS
/

&IONS

&CELL
/

ATOMIC_SPECIES
{self.input.atoms.get_scf_atomic_species()}

CELL_PARAMETERS angstrom
{self.input.atoms.get_scf_cell()}

ATOMIC_POSITIONS angstrom 
{self.input.atoms.get_scf_atomic_positions()}

{self.input.scf.get_kgrid()}
'''
        self.job_scf: str = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.scf.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.scf.job_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.scf.job_desc)} < scf.in &> scf.in.out

cp ./tmp/struct.save/data-file-schema.xml ./scf.xml
'''
    
        self.jobs = [
            'job_scf.sh',
        ]

    def create(self):
        write_str_2_f('scf.in', self.input_scf)
        write_str_2_f('job_scf.sh', self.job_scf)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_scf.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'scf.in',
            'job_scf.sh',
            'tmp',
            'scf.xml',
            '*.pkl',
            '*.xsf',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf scf.in')
        os.system('rm -rf job_scf.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf scf.in.out')
        os.system('rm -rf uc_*.txt')
        os.system('rm -rf sc_*.txt')
        os.system('rm -rf scf.xml')
        os.system('rm -rf bandpath.json')
        os.system('rm -rf kgrid.inp kgrid.log kgrid.out')
        os.system('rm -rf *.xsf')
        os.system('rm -rf ONCVPSP')
#endregion
