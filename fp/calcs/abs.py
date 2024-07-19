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
class Absorption:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_absorption = \
f'''# Q-points
exciton_Q_shift 2 {self.input.absorption.Qshift[0]:15.10f} {self.input.absorption.Qshift[1]:15.10f} {self.input.absorption.Qshift[2]:15.10f}
use_symmetries_coarse_grid
no_symmetries_fine_grid
no_symmetries_shifted_grid

# Bands
number_val_bands_coarse {self.input.absorption.val_bands_coarse}
number_cond_bands_coarse {self.input.absorption.cond_bands_coarse}
number_val_bands_fine {self.input.absorption.val_bands_fine}
number_cond_bands_fine {self.input.absorption.cond_bands_fine}
degeneracy_check_override
#spinor

# Options
diagonalization
#use_velocity
use_momentum
polarization {self.input.absorption.pol_dir[0]:15.10f} {self.input.absorption.pol_dir[1]:15.10f} {self.input.absorption.pol_dir[2]:15.10f}
eqp_co_corrections
dump_bse_hamiltonian

# IO
use_wfn_hdf5

# Output
energy_resolution 0.01
write_eigenvectors {self.input.absorption.num_evec}
'''
        
        self.job_absorption = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.absorption.job_desc)}

ln -sf {self.input.absorption.wfn_co_link} WFN_co.h5 
ln -sf {self.input.absorption.wfnq_co_link} WFNq_co.h5 
ln -sf {self.input.absorption.wfn_fi_link} WFN_fi.h5 
ln -sf {self.input.absorption.wfnq_fi_link} WFNq_fi.h5 
ln -sf eqp1.dat eqp_co.dat 
{self.input.scheduler.get_sched_mpi_prefix(self.input.absorption.job_desc)}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''

        self.jobs = [
            'job_absorption.sh',
        ]

    def create(self):
        write_str_2_f('absorption.inp', self.input_absorption)
        write_str_2_f('job_absorption.sh', self.job_absorption)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_absorption.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'absorption.inp*',
            'eigenvalues.dat',
            'eigenvalues_noeh.dat',
            'absorption_eh.dat',
            'absorption_noeh.dat',
            'eigenvectors.h5',
            'bandstructure_absorption.dat',
            'hbse*.h5',
            'job_absorption.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf absorption.inp')
        os.system('rm -rf job_absorption.sh')
        
        os.system('rm -rf ./WFN_co.h5')
        os.system('rm -rf ./WFNq_co.h5')
        os.system('rm -rf ./WFN_fi.h5')
        os.system('rm -rf ./WFNq_fi.h5')
        os.system('rm -rf eigenvalues.dat')
        os.system('rm -rf eigenvalues_noeh.dat')
        os.system('rm -rf absorption_eh.dat')
        os.system('rm -rf absorption_noeh.dat')
        os.system('rm -rf dvmat_norm.dat')
        os.system('rm -rf dcmat_norm.dat')
        os.system('rm -rf eqp_co.dat')
        os.system('rm -rf eqp.dat')
        os.system('rm -rf eqp_q.dat')
        os.system('rm -rf bandstructure_absorption.dat')
        os.system('rm -rf eigenvectors.h5')
        os.system('rm -rf hbse*.h5')
        os.system('rm -rf x.dat')
        os.system('rm -rf epsdiag.dat')
        os.system('rm -rf dtmat')
        os.system('rm -rf vmtxel')
        os.system('rm -rf absorption.inp.out')

class Plotxct:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.input_plotxct = \
f'''# Cell parameters.
hole_position {self.input.plotxct.get_hole_position_str()}
supercell_size {self.input.plotxct.get_supercell_size_str()}

# Q-points. 
# q_shift
no_symmetries_fine_grid
no_symmetries_shifted_grid

# Bands and state. 
plot_spin 1
plot_state {self.input.plotxct.state}
#spinor
#electron_spin 1
#hole_spin 2

# Output. 

# IO
use_wfn_hdf5
'''
        
        self.job_plotxct = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.plotxct.job_desc)}


ln -sf {self.input.plotxct.wfn_fi_link} WFN_fi.h5 
ln -sf {self.input.plotxct.wfnq_fi_link} WFNq_fi.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.plotxct.job_desc)}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''

        self.jobs = [
            'job_plotxct.sh'
        ]

    def create(self):
        write_str_2_f('plotxct.inp', self.input_plotxct)
        write_str_2_f('job_plotxct.sh', self.job_plotxct)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_plotxct.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'plotxct.inp*',
            'job_plotxct.sh',
            'plotxct_elec.xsf',
            'plotxct_hole.xsf',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf plotxct.inp')
        os.system('rm -rf job_plotxct.sh')
        
        os.system('rm -rf *.a3Dr')
        os.system('rm -rf plotxct.xsf')
        os.system('rm -rf plotxct_elec.xsf')
        os.system('rm -rf plotxct.inp.out')
        os.system('rm -rf WFN_fi.h5')
        os.system('rm -rf WFNq_fi.h5')
#endregion
