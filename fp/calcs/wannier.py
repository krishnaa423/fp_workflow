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
class Wannier:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_wfnwan = \
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
nbnd={self.input.wannier.num_bands}
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

{self.input.wannier.get_kpoints_qe()}
'''
        
        self.job_wfnwan = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wannier.job_wfnwan_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wannier.job_wfnwan_desc)}pw.x < wfnwan.in &> wfnwan.in.out 
'''
        
        self.input_wan = \
f'''# Structure. 
{self.input.wannier.get_unit_cell_cart()}

{self.input.wannier.get_atoms_cart()}

# kpoints. 
{self.input.wannier.get_mpgrid()}

{self.input.wannier.get_kpoints()}

# Bands. 
num_bands = {self.input.wannier.num_bands}
num_wann = {self.input.wannier.num_wann}

# Options. 
auto_projections = .true. 

# Output. 
wannier_plot = .true. 
write_hr = .true.
write_u_matrices = .true. 
'''
        
        self.job_wanpp = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wannier.job_wan_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wannier.job_wan_desc)}wannier90.x {self.input.scheduler.get_sched_mpi_infix(self.input.wannier.job_wan_desc)} -pp wan &> wan.win.pp.out
'''
        
        self.input_pw2wan = \
f'''&INPUTPP
outdir='./tmp'
prefix='struct'
seedname='wan'
write_amn=.true.
write_mmn=.true.
write_unk=.true.
scdm_proj=.true.
/
'''
        
        self.job_pw2wan = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wannier.job_pw2wan_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wannier.job_pw2wan_desc)}pw2wannier90.x < pw2wan.in &> pw2wan.in.out 
'''
        
        self.job_wan = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wannier.job_wan_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wannier.job_wan_desc)}wannier90.x {self.input.scheduler.get_sched_mpi_infix(self.input.wannier.job_wan_desc)} wan  &> wan.win.out 
'''
        
        self.jobs = [
            'job_wfnwan.sh',
            'job_wanpp.sh',
            'job_pw2wan.sh',
            'job_wan.sh',
        ]

    def create(self):
        write_str_2_f('wfnwan.in', self.input_wfnwan)
        write_str_2_f('job_wfnwan.sh', self.job_wfnwan)
        write_str_2_f('wan.win', self.input_wan)
        write_str_2_f('job_wanpp.sh', self.job_wanpp)
        write_str_2_f('pw2wan.in', self.input_pw2wan)
        write_str_2_f('job_pw2wan.sh', self.job_pw2wan)
        write_str_2_f('job_wan.sh', self.job_wan)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_wfnwan.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_wanpp.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_pw2wan.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_wan.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wan*',
            'wfnwan.in*',
            'job_wfnwan.sh',
            'job_wanpp.sh',
            'job_pw2wan.sh',
            'job_wan.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wan*')
        os.system('rm -rf UNK*')
        
        os.system('rm -rf pw2wan.in')
        os.system('rm -rf pw2wan.in.out')
        os.system('rm -rf wfnwan.in')
        os.system('rm -rf wfnwan.in.out')
        
        
        os.system('rm -rf job_wfnwan.sh')
        os.system('rm -rf job_wanpp.sh')
        os.system('rm -rf job_pw2wan.sh')
        os.system('rm -rf job_wan.sh')
#endregion
