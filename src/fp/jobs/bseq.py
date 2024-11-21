#region: Modules.
from fp.inputs.input_main import *
from fp.io.strings import *
from fp.flows.run import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class BseqJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

        self.jobs = [
            'job_bseq.sh',
        ]

        self.bseq_for_xctph_link_str: str = ''

    def get_plotxct_strings(self, Qpt):

        plotxct_spinorbit_extra_args = \
f'''
# Spin-orbit extra args.
spinor 
hole_spin 1
electron_spin 2
'''

        input_plotxct = \
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

{plotxct_spinorbit_extra_args if self.input.scf.is_spinorbit else ""}
{self.input.plotxct.extra_args if self.input.plotxct.extra_args is not None else ""}
'''
        job_plotxct = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.plotxct.job_desc)}

ln -sf {self.input.plotxct.wfn_fi_link} WFN_fi.h5 
ln -sf {self.input.plotxct.wfnq_fi_link} WFNq_fi.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.plotxct.job_desc)}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        
        return input_plotxct, job_plotxct

    def get_kernel_strings(self, Qpt):
        input_ker = \
f'''# Q-points
exciton_Q_shift 2 {Qpt[0]:15.10f} {Qpt[1]:15.10f} {Qpt[2]:15.10f}
use_symmetries_coarse_grid

# Bands 
number_val_bands {self.input.absorption.val_bands_coarse}
number_cond_bands {self.input.absorption.cond_bands_coarse}
#spinor

# Options
#extended_kernel

# IO. 
use_wfn_hdf5
{self.input.kernel.extra_args if self.input.kernel.extra_args is not None else ""}
'''
        
        job_ker = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.kernel.job_desc)}

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../{self.input.absorption.wfn_co_link} WFN_co.h5
ln -sf ../../{self.input.absorption.wfnq_co_link} WFNq_co.h5
{self.input.scheduler.get_sched_mpi_prefix(self.input.kernel.job_desc)}kernel.cplx.x &> kernel.inp.out
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
polarization {self.input.absorption.pol_dir[0]:15.10f} {self.input.absorption.pol_dir[1]:15.10f} {self.input.absorption.pol_dir[2]:15.10f}
eqp_co_corrections
dump_bse_hamiltonian

# IO
use_wfn_hdf5

# Output
energy_resolution 0.01
write_eigenvectors {self.input.absorption.num_evec}
{self.input.absorption.extra_args if self.input.absorption.extra_args is not None else ""}
'''
        
        job_abs = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.absorption.job_desc)}

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../eqp1.dat eqp_co.dat 
ln -sf ../../bsemat.h5 ./
ln -sf ../../{self.input.absorption.wfn_co_link} WFN_co.h5 
ln -sf ../../{self.input.absorption.wfnq_co_link} WFNq_co.h5 
ln -sf ../../{self.input.absorption.wfn_fi_link} WFN_fi.h5 
ln -sf ../../{self.input.absorption.wfnq_fi_link} WFNq_fi.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.absorption.job_desc)}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        
        return input_abs, job_abs

    def create_inputs_bseq(self):
        
        self.bseq_for_xctph_link_str: str = '\n\n\n\n'

        os.system('mkdir -p ./bseq')
        os.system('mkdir -p ./bseq_for_xctph')
        os.chdir('./bseq')

        Qpts = self.input.bseq.get_Qpts()

        for Qpt_idx, Qpt in enumerate(Qpts):
            Qpt0 = f'{Qpt[0]:15.10f}'.strip()
            Qpt1 = f'{Qpt[1]:15.10f}'.strip()
            Qpt2 = f'{Qpt[2]:15.10f}'.strip()
            dir_name = f'Q_{Qpt0}_{Qpt1}_{Qpt2}'
            os.system(f'mkdir -p {dir_name}')
            self.bseq_for_xctph_link_str += f'ln -sf ../bseq/{dir_name} ./bseq_for_xctph/Q_{str(Qpt_idx).strip()}\n'
            os.chdir(f'./{dir_name}')
            
            inp_ker, job_ker = self.get_kernel_strings(Qpt)
            write_str_2_f('kernel.inp', inp_ker)
            write_str_2_f('job_kernel.sh', job_ker)

            inp_abs, job_abs = self.get_absorption_strings(Qpt)
            write_str_2_f('absorption.inp', inp_abs)
            write_str_2_f('job_absorption.sh', job_abs)


            inp_plotxct, job_plotxct = self.get_plotxct_strings(Qpt)
            write_str_2_f('plotxct.inp', inp_plotxct)
            write_str_2_f('job_plotxct.sh', job_plotxct)

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


        kernel_commands = \
f'''ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../{self.input.kernel.wfn_co_link} WFN_co.h5
ln -sf ../../{self.input.kernel.wfnq_co_link} WFNq_co.h5
{self.input.scheduler.get_sched_mpi_prefix(self.input.bseq.job_desc)}kernel.cplx.x &> kernel.inp.out
'''
        
        absorption_commands = \
f'''ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../eqp1.dat eqp_co.dat 
ln -sf ../../bsemat.h5 ./
ln -sf ../../{self.input.absorption.wfn_co_link} WFN_co.h5 
ln -sf ../../{self.input.absorption.wfnq_co_link} WFNq_co.h5 
ln -sf ../../{self.input.absorption.wfn_fi_link} WFN_fi.h5 
ln -sf ../../{self.input.absorption.wfnq_fi_link} WFNq_fi.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.bseq.job_desc)}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        
        plotxct_commands = \
f'''ln -sf ../../{self.input.plotxct.wfn_fi_link} WFN_fi.h5 
ln -sf ../../{self.input.plotxct.wfnq_fi_link} WFNq_fi.h5 
{self.input.scheduler.get_sched_mpi_prefix(self.input.bseq.job_desc)}plotxct.cplx.x &> plotxct.inp.out 
volume.py ../../scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        
        folder_variable = '${folders[$i]}'
        kpt_variable = '${i}'

        # Add the looping block.
        job_bseq += \
f'''
rm -rf ./bseq.out
touch ./bseq.out

{self.bseq_for_xctph_link_str}

for (( i=$start; i<$stop; i++ )); do
    cd {folder_variable}

    echo -e "\\n\\n\\n" >> ../../bseq.out

    echo "Running {kpt_variable} th kpoint" >> ../../bseq.out
    echo "Entering folder {folder_variable}" >> ../../bseq.out
    
    echo "Starting kernel for {folder_variable}" >> ../../bseq.out
{kernel_commands}
    echo "Done kernel for {folder_variable}" >> ../../bseq.out

    echo "Starting absorption for {folder_variable}" >> ../../bseq.out
{absorption_commands}
    echo "Done absorption for {folder_variable}" >> ../../bseq.out

    echo "Starting plotxct for {folder_variable}" >> ../../bseq.out
{plotxct_commands}
    echo "Done plotxct for {folder_variable}" >> ../../bseq.out
    cd ../../

    echo "Exiting folder {folder_variable}" >> ../../bseq.out
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
        os.system('rm -rf ./bseq.out')
        os.system('rm -rf ./bseq_for_xctph')
        os.system('rm -rf ./job_bseq.sh')
#endregion
