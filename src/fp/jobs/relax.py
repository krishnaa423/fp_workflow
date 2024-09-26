#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import *
import os 
from fp.flows.run import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class RelaxJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.save_final_cell_parameters_str = \
"awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /CELL_PARAMETERS/ {cell_flag=1; next} /End final coordinates/ {end_flag=0} end_flag && cell_flag {print; if (length==0) cell_flag=0 }' relax.in.out > relaxed_cell_parameters.txt"
        
        # Skipping the first field. 
        # "awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /ATOMIC_POSITIONS/ {pos_flag=1} /End final coordinates/ {end_flag=0}  end_flag && pos_flag {print}' relax.in.out > relaxed_atomic_positions.txt"
        self.save_final_atomic_positions_str = \
"awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /ATOMIC_POSITIONS/ {pos_flag=1; next} /End final coordinates/ {end_flag=0}  end_flag && pos_flag { print $2, $3, $4 }' relax.in.out > relaxed_atomic_positions.txt"

        self.input_relax: str = \
f'''&CONTROL
outdir='./tmp'
prefix='struct'
pseudo_dir='./ONCVPSP/sg15'
calculation='{self.input.relax.calc_str()}'
tprnfor=.true.
{self.input.relax.extra_control_args if self.input.relax.extra_control_args is not None else ""}
/

&SYSTEM
ibrav=0
occupations='from_input'
ntyp={self.input.atoms.get_ntyp()}
nat={self.input.atoms.get_nat()}
nbnd={self.input.relax.get_nbnd()}
ecutwfc={self.input.scf.ecutwfc}
{"" if self.input.scf.is_spinorbit else "!"}noncolin=.true.
{"" if self.input.scf.is_spinorbit else "!"}lspinorb=.true. 
{self.input.relax.extra_system_args if self.input.relax.extra_system_args is not None else ""}
/

&ELECTRONS
{self.input.relax.extra_electrons_args if self.input.relax.extra_electrons_args is not None else ""}
/

&IONS
/

&CELL
/

{self.input.relax.get_occupations_str()}

ATOMIC_SPECIES
{self.input.atoms.get_qe_scf_atomic_species(self.input.scf.is_spinorbit)}

CELL_PARAMETERS angstrom
{self.input.atoms.get_scf_cell()}

ATOMIC_POSITIONS angstrom 
{self.input.atoms.get_qe_scf_atomic_positions()}

{self.input.scf.get_kgrid()}
'''
        self.job_relax: str = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.scf.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.relax.job_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.relax.job_desc)} < relax.in &> relax.in.out

cp ./tmp/struct.save/data-file-schema.xml ./relax.xml

# Copy the end atomic positions and cell parameters (if vc-relax).
{self.save_final_cell_parameters_str}
{self.save_final_atomic_positions_str}
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
