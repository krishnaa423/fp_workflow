#region: Modules.
from fp.inputs import *
from fp.io.strings import *
from fp.flows.run import *
import re 
import subprocess
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Dryrun:
    def __init__(
        self,
        atoms: AtomsInput,
    ):
        self.atoms: AtomsInput = atoms

        self.input_dryrun: str = \
f'''&CONTROL
outdir='./tmp'
prefix='struct'
pseudo_dir='./ONCVPSP/sg15'
calculation='md'
nstep=0
tprnfor=.true.
/

&SYSTEM
ibrav=0
ntyp={self.atoms.get_ntyp()}
nat={self.atoms.get_nat()}
ecutwfc=20.0
!noncolin=.true.
!lspinorb=.true.
/

&ELECTRONS
/

&IONS

&CELL
/

ATOMIC_SPECIES
{self.atoms.get_scf_atomic_species()}

CELL_PARAMETERS angstrom
{self.atoms.get_scf_cell()}

ATOMIC_POSITIONS angstrom 
{self.atoms.get_scf_atomic_positions()}

K_POINTS automatic 
2 2 2 0 0 0
'''
        # The dryrun jobs will not be in parallel. Just to get some info.
        self.job_dryrun: str = \
f'''#!/bin/bash

pw.x < dryrun.in &> dryrun.in.out

cp ./tmp/struct.save/data-file-schema.xml ./dryrun.xml
'''
    
    def create(self):
        write_str_2_f('dryrun.in', self.input_dryrun)
        write_str_2_f('job_dryrun.sh', self.job_dryrun)

    def run(self, total_time):
        subprocess.run('./job_dryrun.sh')

        return 0.0

    def save(self, folder):
        pass 

    def remove(self):
        os.system('rm -rf dryrun.in')
        os.system('rm -rf job_dryrun.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf dryrun.in.out')
        os.system('rm -rf dryrun.xml')

    def get_max_val(self):
        with open('dryrun.in.out', 'r') as r: txt = r.read()
        pattern  = r'number of Kohn-Sham states=(.*)\n'
        num_kohn_sham_states = int(re.findall(pattern, txt)[0])
        
        return num_kohn_sham_states
#endregion
