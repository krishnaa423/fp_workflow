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
class Epw:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_epw = \
f'''&INPUTEPW
outdir='./tmp'
prefix='struct'

! kpoints.
nk1={self.input.epw.kgrid_coarse[0]}
nk2={self.input.epw.kgrid_coarse[1]}
nk3={self.input.epw.kgrid_coarse[2]}
nq1={self.input.epw.qgrid_coarse[0]}
nq2={self.input.epw.qgrid_coarse[1]}
nq3={self.input.epw.qgrid_coarse[2]}
nkf1={self.input.epw.kgrid_fine[0]}
nkf2={self.input.epw.kgrid_fine[1]}
nkf3={self.input.epw.kgrid_fine[2]}
nqf1={self.input.epw.qgrid_fine[0]}
nqf2={self.input.epw.qgrid_fine[1]}
nqf3={self.input.epw.qgrid_fine[2]}

! Bands. 
nbndsub={self.input.epw.bands}
{self.input.epw.get_skipped_bands_str()}

! elph. 
dvscf_dir='./save' 
elph=.true. 
epbwrite=.true. 
epbread=.false.
!prtgkk=.true.

! wannier. 
wannierize=.true. 
proj(1)='random'

! others. 
!temps=300.0
!iverbosity=1
/
'''
        
        self.job_epw = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.epw.job_desc)}

{self.input.scheduler.get_sched_mpi_prefix(self.input.epw.job_desc)}{self.input.epw.exec_loc} {self.input.scheduler.get_sched_mpi_infix(self.input.epw.job_desc)} < epw.in &> epw.in.out 
cp ./tmp/struct_elph* ./
'''
        
        self.jobs = [
            'job_epw.sh',
        ]

    def create(self):
        write_str_2_f('epw.in', self.input_epw)
        write_str_2_f('job_epw.sh', self.job_epw)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_epw.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'epw.in*',
            'struct_elph*',
            'job_epw.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf epw.in')
        os.system('rm -rf job_epw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf ./struct*')
        os.system('rm -rf ./decay*')
        os.system('rm -rf ./struct_elph*')
        os.system('rm -rf EPW.bib')
        os.system('rm -rf epwdata.fmt')
        os.system('rm -rf selecq.fmt')
        os.system('rm -rf vmedata.fmt')
        os.system('rm -rf crystal.fmt')
        os.system('rm -rf epw.in.out')
#endregion
