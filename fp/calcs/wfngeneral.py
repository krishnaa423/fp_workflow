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
wfng_file='WFN'
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
cp ./tmp/WFN ./
cp ./tmp/RHO ./
cp ./tmp/VXC ./
cp ./tmp/VSC ./
cp ./tmp/VKB ./
'''

        self.input_parabands = \
f'''input_wfn_file WFN
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
        pass 

    def remove(self):
        os.system('rm -rf wfn.in')
        os.system('rm -rf job_wfn.sh')
        os.system('rm -rf wfn_pw2bgw.in')
        os.system('rm -rf job_wfn_pw2bgw.sh')
        os.system('rm -rf parabands.inp')
        os.system('rm -rf job_parabands.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf ./WFN')
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
wfng_file='WFNq'
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
cp ./tmp/WFNq ./
wfn2hdf.x BIN WFNq WFNq.h5 
'''

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
        pass 

    def remove(self):
        os.system('rm -rf wfnq.in')
        os.system('rm -rf job_wfnq.sh')
        os.system('rm -rf wfnq_pw2bgw.in')
        os.system('rm -rf job_wfnq_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFNq')
        os.system('rm -rf WFNq.h5')
        os.system('rm -rf wfnq.in.out')
        os.system('rm -rf wfnq_pw2bgw.in.out')


#endregion
