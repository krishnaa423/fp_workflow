#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import *
import os 
from fp.flows import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Relax:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_relax: str = \
f'''&CONTROL
outdir='./tmp'
prefix='struct'
pseudo_dir='./ONCVPSP/sg15'
calculation='{self.input.relax.calc_str()}'
tprnfor=.true.
/

&SYSTEM
ibrav=0
occupations='from_input'
ntyp={self.input.atoms.get_ntyp()}
nat={self.input.atoms.get_nat()}
nbnd={self.input.relax.get_nbnd()}
ecutwfc={self.input.scf.ecutwfc}
!noncolin=.true.
!lspinorb=.true. 
/

&ELECTRONS
/

&IONS

&CELL
/

{self.input.relax.get_occupations_str()}

ATOMIC_SPECIES
{self.input.atoms.get_scf_atomic_species()}

CELL_PARAMETERS angstrom
{self.input.atoms.get_scf_cell()}

ATOMIC_POSITIONS angstrom 
{self.input.atoms.get_scf_atomic_positions()}

{self.input.scf.get_kgrid()}
'''
        self.job_relax: str = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.scf.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.relax.job_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.relax.job_desc)} < relax.in &> relax.in.out

cp ./tmp/struct.save/data-file-schema.xml ./relax.xml
'''
    
        self.jobs = [
            'job_relax.sh',
        ]

    def create(self):
        write_str_2_f('relax.in', self.input_relax)
        write_str_2_f('job_relax.sh', self.job_relax)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_relax.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'relax.in',
            'job_relax.sh',
            'relax.in.out',
            'relax.xml',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf relax.in')
        os.system('rm -rf job_relax.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf relax.in.out')
        os.system('rm -rf relax.xml')
#endregion
