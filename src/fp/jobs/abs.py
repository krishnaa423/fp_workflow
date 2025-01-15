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
class AbsorptionJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['abs']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['abs']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['abs']['job_info'])

    def set_inputs_str(self):
        # Base.
        input_absorption_dict: dict = {
            'maps': {
                'exciton_Q_shift': [
                    2,
                    self.input_dict['abs']['Qshift'][0],
                    self.input_dict['abs']['Qshift'][1],
                    self.input_dict['abs']['Qshift'][2],
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
        self.input_absorption: str = BgwInputFile.write_general(input_absorption_dict)

    def set_jobs_str(self):
        self.job_absorption = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

ln -sf {self.input_dict['abs']['wfnco_link']} WFN_co.h5 
ln -sf {self.input_dict['abs']['wfnqco_link']} WFNq_co.h5 
ln -sf {self.input_dict['abs']['wfnfi_link']} WFN_fi.h5 
ln -sf {self.input_dict['abs']['wfnqfi_link']} WFNq_fi.h5 
ln -sf eqp1.dat eqp_co.dat 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''

        self.jobs = [
            './job_absorption.sh',
        ]

    def create(self):
        write_str_2_f('absorption.inp', self.input_absorption)
        write_str_2_f('job_absorption.sh', self.job_absorption)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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

class PlotxctJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['plotxct']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['plotxct']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['plotxct']['job_info'])

    def set_inputs_str(self):
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
        self.input_plotxct: str = BgwInputFile.write_general(input_plotxct_dict)

    def set_jobs_str(self):
        self.job_plotxct = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}


ln -sf {self.input_dict['plotxct']['wfnfi_link']} WFN_fi.h5 
ln -sf {self.input_dict['plotxct']['wfnqfi_link']} WFNq_fi.h5 
{self.scheduler.get_sched_mpi_prefix(self.job_info)}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''

        self.jobs = [
            './job_plotxct.sh'
        ]

    def create(self):
        write_str_2_f('plotxct.inp', self.input_plotxct)
        write_str_2_f('job_plotxct.sh', self.job_plotxct)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
