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
class Wfn:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_wfn = \
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
nbnd={self.input.wfn.bands}
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

{self.input.wfn.get_kgrid_dft_string()}
'''
        
        self.job_wfn = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfn.job_wfn_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfn.job_wfn_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.wfn.job_wfn_desc)} < wfn.in &> wfn.in.out 
'''
        
        self.input_wfn_pw2bgw = \
f'''&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFN_coo'
wfng_kgrid=.true.
wfng_nk1={self.input.wfn.kdim[0]}
wfng_nk2={self.input.wfn.kdim[1]}
wfng_nk3={self.input.wfn.kdim[2]}
wfng_dk1={self.input.wfn.qshift[0]*self.input.wfn.kdim[0]}
wfng_dk2={self.input.wfn.qshift[1]*self.input.wfn.kdim[1]}
wfng_dk3={self.input.wfn.qshift[2]*self.input.wfn.kdim[2]}
rhog_flag=.true.
rhog_file='RHO'
vxcg_flag=.true.
vxcg_file='VXC'
vscg_flag=.true.
vscg_file='VSC'
vkbg_flag=.true.
vkbg_file='VKB'
/
'''

        self.job_wfn_pw2bgw = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfn.job_pw2bgw_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfn.job_pw2bgw_desc)}pw2bgw.x -pd .true. < wfn_pw2bgw.in &> wfn_pw2bgw.in.out 
cp ./tmp/WFN_coo ./
cp ./tmp/RHO ./
cp ./tmp/VXC ./
cp ./tmp/VSC ./
cp ./tmp/VKB ./
'''

        self.input_parabands = \
f'''input_wfn_file WFN_coo
output_wfn_file WFN_parabands.h5 

vsc_file VSC 
vkb_file VKB 

number_bands {self.input.wfn.parabands_bands}
'''
        
        self.job_parabands = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfn.job_parabands_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfn.job_parabands_desc)}parabands.cplx.x &> parabands.inp.out 
'''

        self.jobs = [
            'job_wfn.sh',
            'job_wfn_pw2bgw.sh',
            'job_parabands.sh',
        ]

    def create(self):
        write_str_2_f(f'wfn.in', self.input_wfn)
        write_str_2_f(f'job_wfn.sh', self.job_wfn)
        write_str_2_f(f'wfn_pw2bgw.in', self.input_wfn_pw2bgw)
        write_str_2_f(f'job_wfn_pw2bgw.sh', self.job_wfn_pw2bgw)
        write_str_2_f(f'parabands.inp', self.input_parabands)
        write_str_2_f(f'job_parabands.sh', self.job_parabands)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_wfn.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_wfn_pw2bgw.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_parabands.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wfn.in*',
            'WFN_coo*',
            'WFN_parabands.h5',
            'job_wfn*',
            'parabands.inp*',
            'job_parabands.sh',
            'RHO',
            'VXC',
            'VSC',
            'VKB',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wfn.in')
        os.system('rm -rf job_wfn.sh')
        os.system('rm -rf wfn_pw2bgw.in')
        os.system('rm -rf job_wfn_pw2bgw.sh')
        os.system('rm -rf parabands.inp')
        os.system('rm -rf job_parabands.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf ./WFN_coo')
        os.system('rm -rf ./WFN_parabands.h5')
        os.system('rm -rf ./RHO')
        os.system('rm -rf ./VXC')
        os.system('rm -rf ./VSC')
        os.system('rm -rf ./VKB')
        os.system('rm -rf wfn.in.out')
        os.system('rm -rf wfn_pw2bgw.in.out')
        os.system('rm -rf parabands.inp.out')

class Wfnq:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_wfnq = \
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
nbnd={self.input.wfnq.bands}
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

{self.input.wfnq.get_kgrid_dft_string()}
'''
        
        self.job_wfnq = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfnq.job_wfn_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfnq.job_wfn_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.wfnq.job_wfn_desc)} < wfnq.in &> wfnq.in.out 
'''
        
        self.input_wfnq_pw2bgw = \
f'''&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFNq_coo'
wfng_kgrid=.true.
wfng_nk1={self.input.wfnq.kdim[0]}
wfng_nk2={self.input.wfnq.kdim[1]}
wfng_nk3={self.input.wfnq.kdim[2]}
wfng_dk1={self.input.wfnq.qshift[0]*self.input.wfnq.kdim[0]}
wfng_dk2={self.input.wfnq.qshift[1]*self.input.wfnq.kdim[1]}
wfng_dk3={self.input.wfnq.qshift[2]*self.input.wfnq.kdim[2]}
/
'''

        self.job_wfnq_pw2bgw = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfnq.job_pw2bgw_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfnq.job_pw2bgw_desc)}pw2bgw.x -pd .true. < wfnq_pw2bgw.in &> wfnq_pw2bgw.in.out 
cp ./tmp/WFNq_coo ./
wfn2hdf.x BIN WFNq_coo WFNq_coo.h5 
'''

        self.jobs = [
            'job_wfnq.sh',
            'job_wfnq_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f(f'wfnq.in', self.input_wfnq)
        write_str_2_f(f'job_wfnq.sh', self.job_wfnq)
        write_str_2_f(f'wfnq_pw2bgw.in', self.input_wfnq_pw2bgw)
        write_str_2_f(f'job_wfnq_pw2bgw.sh', self.job_wfnq_pw2bgw)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_wfnq.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_wfnq_pw2bgw.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wfnq.in*',
            'WFNq_coo*',
            'job_wfnq*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wfnq.in')
        os.system('rm -rf job_wfnq.sh')
        os.system('rm -rf wfnq_pw2bgw.in')
        os.system('rm -rf job_wfnq_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFNq_coo')
        os.system('rm -rf WFNq_coo.h5')
        os.system('rm -rf wfnq.in.out')
        os.system('rm -rf wfnq_pw2bgw.in.out')

class Wfnfi:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_wfnfi = \
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
nbnd={self.input.wfnfi.bands}
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

{self.input.wfnfi.get_kgrid_dft_string()}
'''
        
        self.job_wfnfi = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfnfi.job_wfn_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfnfi.job_wfn_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.wfnfi.job_wfn_desc)} < wfnfi.in &> wfnfi.in.out 
'''
        
        self.input_wfnfi_pw2bgw = \
f'''&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFN_fii'
wfng_kgrid=.true.
wfng_nk1={self.input.wfnfi.kdim[0]}
wfng_nk2={self.input.wfnfi.kdim[1]}
wfng_nk3={self.input.wfnfi.kdim[2]}
wfng_dk1={self.input.wfnfi.qshift[0]*self.input.wfnfi.kdim[0]}
wfng_dk2={self.input.wfnfi.qshift[1]*self.input.wfnfi.kdim[1]}
wfng_dk3={self.input.wfnfi.qshift[2]*self.input.wfnfi.kdim[2]}
/
'''

        self.job_wfnfi_pw2bgw = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfnfi.job_pw2bgw_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfnfi.job_pw2bgw_desc)}pw2bgw.x -pd .true. < wfnfi_pw2bgw.in &> wfnfi_pw2bgw.in.out 
cp ./tmp/WFN_fii ./
wfn2hdf.x BIN WFN_fii WFN_fii.h5 
'''

        self.jobs = [
            'job_wfnfi.sh',
            'job_wfnfi_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f(f'wfnfi.in', self.input_wfnfi)
        write_str_2_f(f'job_wfnfi.sh', self.job_wfnfi)
        write_str_2_f(f'wfnfi_pw2bgw.in', self.input_wfnfi_pw2bgw)
        write_str_2_f(f'job_wfnfi_pw2bgw.sh', self.job_wfnfi_pw2bgw)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_wfnfi.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_wfnfi_pw2bgw.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wfnfi.in*',
            'WFN_fii*',
            'job_wfnfi*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')


    def remove(self):
        os.system('rm -rf wfnfi.in')
        os.system('rm -rf job_wfnfi.sh')
        os.system('rm -rf wfnfi_pw2bgw.in')
        os.system('rm -rf job_wfnfi_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFN_fii')
        os.system('rm -rf WFN_fii.h5')
        os.system('rm -rf wfnfi.in.out')
        os.system('rm -rf wfnfi_pw2bgw.in.out')

class Wfnqfi:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_wfnqfi = \
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
nbnd={self.input.wfnqfi.bands}
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

{self.input.wfnqfi.get_kgrid_dft_string()}
'''
        
        self.job_wfnqfi = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfnqfi.job_wfn_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfnqfi.job_wfn_desc)}pw.x {self.input.scheduler.get_sched_mpi_infix(self.input.wfnqfi.job_wfn_desc)} < wfnqfi.in &> wfnqfi.in.out 
'''
        
        self.input_wfnqfi_pw2bgw = \
f'''&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFNq_fii'
wfng_kgrid=.true.
wfng_nk1={self.input.wfnqfi.kdim[0]}
wfng_nk2={self.input.wfnqfi.kdim[1]}
wfng_nk3={self.input.wfnqfi.kdim[2]}
wfng_dk1={self.input.wfnqfi.qshift[0]*self.input.wfnqfi.kdim[0]}
wfng_dk2={self.input.wfnqfi.qshift[1]*self.input.wfnqfi.kdim[1]}
wfng_dk3={self.input.wfnqfi.qshift[2]*self.input.wfnqfi.kdim[2]}
/
'''

        self.job_wfnqfi_pw2bgw = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.wfnqfi.job_pw2bgw_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.wfnqfi.job_pw2bgw_desc)}pw2bgw.x -pd .true. < wfnqfi_pw2bgw.in &> wfnqfi_pw2bgw.in.out 
cp ./tmp/WFNq_fii ./
wfn2hdf.x BIN WFNq_fii WFNq_fii.h5 
'''

        self.jobs = [
            'job_wfnqfi.sh',
            'job_wfnqfi_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f(f'wfnqfi.in', self.input_wfnqfi)
        write_str_2_f(f'job_wfnqfi.sh', self.job_wfnqfi)
        write_str_2_f(f'wfnqfi_pw2bgw.in', self.input_wfnqfi_pw2bgw)
        write_str_2_f(f'job_wfnqfi_pw2bgw.sh', self.job_wfnqfi_pw2bgw)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_wfnqfi.sh', self.input, total_time)
        total_time = run_and_wait_command('./job_wfnqfi_pw2bgw.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'wfnqfi.in*',
            'WFNq_fii*',
            'job_wfnqfi*',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf wfnqfi.in')
        os.system('rm -rf job_wfnqfi.sh')
        os.system('rm -rf wfnqfi_pw2bgw.in')
        os.system('rm -rf job_wfnqfi_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFNq_fii')
        os.system('rm -rf WFNq_fii.h5')
        os.system('rm -rf wfnqfi.in.out')
        os.system('rm -rf wfnqfi_pw2bgw.in.out')

#endregion
