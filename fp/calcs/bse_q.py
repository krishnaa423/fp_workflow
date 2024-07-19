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
class Bseq:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.jobs = [
            'job_bseq.sh',
        ]

    def get_kernel_strings(self, Qpt):
        input_ker = \
f'''# Q-points
exciton_Q_shift 2 {Qpt[0]:15.10f} {Qpt[1]:15.10f} {Qpt[2]:15.10f}
use_symmetries_coarse_grid

# Bands 
number_val_bands {self.input.bseq.val_bands_coarse}
number_cond_bands {self.input.bseq.cond_bands_coarse}
#spinor

# Options
#extended_kernel

# IO. 
use_wfn_hdf5
'''
        
        job_ker = \
f'''#!/bin/bash

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../{self.input.bseq.wfn_co_link} WFN_co.h5
ln -sf ../../{self.input.bseq.wfnq_co_link} WFNq_co.h5
{self.input.scheduler.get_sched_mpi_prefix(self.input.bseq.job_desc)}kernel.cplx.x &> kernel.inp.out
'''
        
        return input_ker, job_ker

    def get_absorption_strings(self, Qpt):
        input_abs = \
f'''# Q-points
exciton_Q_shift 2 {Qpt[0]:15.10f} {Qpt[1]:15.10f} {Qpt[2]:15.10f}
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
polarization {self.input.bseq.pol_dir[0]:15.10f} {self.input.bseq.pol_dir[1]:15.10f} {self.input.bseq.pol_dir[2]:15.10f}
eqp_co_corrections
dump_bse_hamiltonian

# IO
use_wfn_hdf5

# Output
energy_resolution 0.01
write_eigenvectors {self.input.bseq.num_evec}
'''
        
        job_abs = \
f'''#!/bin/bash

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../eqp1.dat eqp_co.dat 
ln -sf ../../bsemat.h5 ./
ln -sf ../../{self.input.bseq.wfn_co_link} WFN_co.h5 
ln -sf ../../{self.input.bseq.wfnq_co_link} WFNq_co.h5 
ln -sf ../../{self.input.bseq.wfn_fi_link} WFN_fi.h5 
ln -sf ../../{self.input.bseq.wfnq_fi_link} WFNq_fi.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.bseq.job_desc)}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        
        return input_abs, job_abs

    def create_inputs_bseq(self):
        
        os.system('mkdir -p ./bseq')
        os.chdir('./bseq')

        Qpts = self.input.bseq.get_Qpts()

        for Qpt in Qpts:
            Qpt0 = f'{Qpt[0]:15.10f}'.strip()
            Qpt1 = f'{Qpt[1]:15.10f}'.strip()
            Qpt2 = f'{Qpt[2]:15.10f}'.strip()
            dir_name = f'Q_{Qpt0}_{Qpt1}_{Qpt2}'
            os.system(f'mkdir -p {dir_name}')
            os.chdir(f'./{dir_name}')
            
            inp_ker, job_ker = self.get_kernel_strings(Qpt)
            write_str_2_f('kernel.inp', inp_ker)
            write_str_2_f('job_kernel.sh', job_ker)

            inp_abs, job_abs = self.get_absorption_strings(Qpt)
            write_str_2_f('absorption.inp', inp_abs)
            write_str_2_f('job_absorption.sh', job_abs)

            os.chdir('../')

        os.chdir('../')

    def create_job_bseq(self):
        '''
        Idea is to create a list with start and stop indices to control execution.
        '''
        Qpts = self.input.bseq.get_Qpts()
        job_bseq = '#!/bin/bash\n'
        job_bseq += f'{self.input.scheduler.get_sched_header(self.input.bseq.job_desc)}\n'

        job_bseq += "start=0\n"
        job_bseq += f"stop={Qpts.shape[0]}\n\n"
        job_bseq += f"size={Qpts.shape[0]}\n\n"

        # Create the list.
        job_bseq += 'folders=('
        for Qpt in Qpts:
            Qpt0 = f'{Qpt[0]:15.10f}'.strip()
            Qpt1 = f'{Qpt[1]:15.10f}'.strip()
            Qpt2 = f'{Qpt[2]:15.10f}'.strip()
            subdir_name = f'Q_{Qpt0}_{Qpt1}_{Qpt2}'
            dir_name = f'"./bseq/{subdir_name}" '
            job_bseq += dir_name
        job_bseq += ')\n\n'

        # Add the looping block.
        job_bseq += \
'''
for (( i=$start; i<$stop; i++ )); do
    cd ${folders[$i]}
    
    echo "Starting kernel for ${folders[$i]}"
    ./job_kernel.sh
    echo "Done kernel for ${folders[$i]}"

    echo "Starting absorption for ${folders[$i]}"
    ./job_absorption.sh
    echo "Done absorption for ${folders[$i]}"
    cd ../../
done
'''

        write_str_2_f('job_bseq.sh', job_bseq)

    def create(self):
        self.create_inputs_bseq()
        self.create_job_bseq()

    def run(self, total_time):
        total_time = run_and_wait_command('./job_bseq.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'bseq',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf ./bseq')
        os.system('rm -rf ./job_bseq.sh')
#endregion
