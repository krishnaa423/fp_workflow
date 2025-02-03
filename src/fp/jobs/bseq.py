#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.bgw import BgwInputFile
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
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.set_job_info()

    def set_job_info(self):
        if isinstance(self.input_dict['bseq']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['bseq']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['bseq']['job_info'])

        self.jobs = [
            './job_bseq.sh',
        ]

    def get_plotxct_strings(self, Qpt):
        # Base.
        input_plotxct_dict: dict = {
            'maps': {
                'hole_position': self.input_dict['plotxct']['hole_position'],
                'supercell_size': self.input_dict['plotxct']['supercell_size'],
                'use_symmetries_fine_grid': '',
                'use_symmetries_shifted_grid': '',
                'plot_spin': 1,
                'plot_state': self.input_dict['plotxct']['xct_state'],
                'use_wfn_hdf5': '',
            },
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_plotxct_dict['maps']['spinor'] = ''
            input_plotxct_dict['maps']['hole_spin'] = 1
            input_plotxct_dict['maps']['electron_spin'] = 2
        args_dict = self.input_dict['plotxct']['args']
        args_type = self.input_dict['plotxct']['args_type']
        input_plotxct_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_plotxct_dict
        )

        # Write.
        input_plotxct: str = BgwInputFile.write_general(input_plotxct_dict)

        job_plotxct = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf ../../{self.input_dict['plotxct']['wfnfi_link']} WFN_fi.h5 
ln -sf ../../{self.input_dict['plotxct']['wfnqfi_link']} WFNq_fi.h5 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        
        return input_plotxct, job_plotxct

    def get_kernel_strings(self, Qpt):
        # Base.
        input_kernel_dict: dict = {
            'maps': {
                'exciton_Q_shift': [
                    2,
                    Qpt[0],
                    Qpt[1],
                    Qpt[2],
                ],
                'use_symmetries_coarse_grid': '',
                'number_val_bands': self.input_dict['abs']['num_val_bands'],
                'number_cond_bands': self.input_dict['abs']['num_cond_bands'],
                'use_wfn_hdf5': '',
            },
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_kernel_dict['maps']['spinor'] = ''
        args_dict = self.input_dict['ker']['args']
        args_type = self.input_dict['ker']['args_type']
        input_kernel_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_kernel_dict
        )

        # Write.
        input_kernel: str = BgwInputFile.write_general(input_kernel_dict)
        
        job_kernel = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../{self.input_dict['abs']['wfnco_link']} WFN_co.h5
ln -sf ../../{self.input_dict['abs']['wfnqco_link']} WFNq_co.h5
{self.scheduler.get_sched_mpi_prefix(self.job_info)}kernel.cplx.x &> kernel.inp.out
'''
        
        return input_kernel, job_kernel

    def get_absorption_strings(self, Qpt):
        # Base.
        input_absorption_dict: dict = {
            'maps': {
                'exciton_Q_shift': [
                    2,
                    Qpt[0],
                    Qpt[1],
                    Qpt[2],
                ],
                'use_symmetries_coarse_grid': '',
                'use_symmetries_fine_grid': '',
                'use_symmetries_shifted_grid': '',
                'number_val_bands_coarse': self.input_dict['abs']['num_val_bands'],
                'number_val_bands_fine': self.input_dict['abs']['num_val_bands'] - 1,
                'number_cond_bands_coarse': self.input_dict['abs']['num_cond_bands'],
                'number_cond_bands_fine': self.input_dict['abs']['num_cond_bands'],
                'degeneracy_check_override': '',
                'diagonalization': '',
                # 'use_elpa': '',  
                'use_momentum': '',  
                # 'use_velocity': '',  
                'polarization': self.input_dict['abs']['pol_dir'],
                'eqp_co_corrections': '',
                'dump_bse_hamiltonian': '',
                'use_wfn_hdf5': '',
                'energy_resolution': 0.1,
                'write_eigenvectors': self.input_dict['abs']['num_evecs'],
            },
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_absorption_dict['maps']['spinor'] = ''
        args_dict = self.input_dict['abs']['args']
        args_type = self.input_dict['abs']['args_type']
        input_absorption_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_absorption_dict
        )

        # Write.
        input_absorption: str = BgwInputFile.write_general(input_absorption_dict)
        
        job_absorption = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../eqp1.dat eqp_co.dat 
#ln -sf ../../bsemat.h5 ./
ln -sf ../../{self.input_dict['abs']['wfnco_link']} WFN_co.h5 
ln -sf ../../{self.input_dict['abs']['wfnqco_link']} WFNq_co.h5 
ln -sf ../../{self.input_dict['abs']['wfnfi_link']} WFN_fi.h5 
ln -sf ../../{self.input_dict['abs']['wfnqfi_link']} WFNq_fi.h5 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        
        return input_absorption, job_absorption

    def create_inputs_bseq(self):
        
        self.bseq_for_xctph_link_str: str = '\n\n\n\n'

        os.system('mkdir -p ./bseq')
        os.system('mkdir -p ./bseq_for_xctph')
        os.chdir('./bseq')

        Qpts = self.input.bseq.get_Qpts(self.input.atoms)

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
        Qpts = self.input.bseq.get_Qpts(self.input.atoms)
        job_bseq = '#!/bin/bash\n'
        job_bseq += f'{self.scheduler.get_sched_header(self.job_info)}\n'

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
f'''    ln -sf ../../epsmat.h5 ./
    ln -sf ../../eps0mat.h5 ./
    ln -sf ../../{self.input_dict['abs']['wfnco_link']} WFN_co.h5
    ln -sf ../../{self.input_dict['abs']['wfnqco_link']} WFNq_co.h5
    {self.scheduler.get_sched_mpi_prefix(self.job_info)}kernel.cplx.x &> kernel.inp.out
'''
        
        absorption_commands = \
f'''    ln -sf ../../epsmat.h5 ./
    ln -sf ../../eps0mat.h5 ./
    ln -sf ../../eqp1.dat eqp_co.dat 
    ln -sf ../../{self.input_dict['abs']['wfnco_link']} WFN_co.h5 
    ln -sf ../../{self.input_dict['abs']['wfnqco_link']} WFNq_co.h5 
    ln -sf ../../{self.input_dict['abs']['wfnfi_link']} WFN_fi.h5 
    ln -sf ../../{self.input_dict['abs']['wfnqfi_link']} WFNq_fi.h5 
    {self.scheduler.get_sched_mpi_prefix(self.job_info)}absorption.cplx.x &> absorption.inp.out
    mv bandstructure.dat bandstructure_absorption.dat
'''
        
        plotxct_commands = \
f'''    ln -sf ../../{self.input_dict['plotxct']['wfnfi_link']} WFN_fi.h5 
    ln -sf ../../{self.input_dict['plotxct']['wfnqfi_link']} WFNq_fi.h5 
    {self.scheduler.get_sched_mpi_prefix(self.job_info)}plotxct.cplx.x &> plotxct.inp.out 
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

LOG_FILE="$(pwd)/bseq.out"
exec &> "$LOG_FILE"

{self.bseq_for_xctph_link_str}

for (( i=$start; i<$stop; i++ )); do
    cd {folder_variable}

    echo -e "\\n\\n\\n"

    echo "Running {kpt_variable} th kpoint"
    echo "Entering folder {folder_variable}"
    
    echo "Starting kernel for {folder_variable}"
{kernel_commands}
    echo "Done kernel for {folder_variable}"

    echo "Starting absorption for {folder_variable}"
{absorption_commands}
    echo "Done absorption for {folder_variable}"

    echo "Starting plotxct for {folder_variable}"
{plotxct_commands}
    echo "Done plotxct for {folder_variable}"
    cd ../../

    echo "Exiting folder {folder_variable}"
done
'''

        write_str_2_f('job_bseq.sh', job_bseq)

    def create(self):
        self.create_inputs_bseq()
        self.create_job_bseq()

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
