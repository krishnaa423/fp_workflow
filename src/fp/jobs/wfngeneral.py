#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
from fp.inputs.bgw import BgwInputFile
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class WfnJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.job_pw2bgw_info: JobProcDesc = None
        self.job_parabands_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['wfn']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['wfn']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['wfn']['job_info'])

        # pw2bgw.
        if isinstance(self.input_dict['wfn']['job_pw2bgw_info'], str):
            self.job_pw2bgw_info = JobProcDesc.from_job_id(
                self.input_dict['wfn']['job_pw2bgw_info'],
                self.input_dict,
            )
        else:
            self.job_pw2bgw_info = JobProcDesc(**self.input_dict['wfn']['job_pw2bgw_info'])

        # parabands.
        if isinstance(self.input_dict['wfn']['job_parabands_info'], str):
            self.job_parabands_info = JobProcDesc.from_job_id(
                self.input_dict['wfn']['job_parabands_info'],
                self.input_dict,
            )
        else:
            self.job_parabands_info = JobProcDesc(**self.input_dict['wfn']['job_parabands_info'])

    def set_inputs_str(self):
        #Base. 
        input_wfn_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'bands',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'nbnd': self.input_dict['wfn']['num_cond_bands'] + self.input_dict['total_valence_bands'],
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                },
                'electrons': {},
                'ions': {},
                'cell': {},
            },
            'blocks': {
                'atomic_species': self.input.atoms.get_qe_scf_atomic_species(),
                'cell_parameters': self.input.atoms.get_qe_scf_cell(),
                'atomic_positions': self.input.atoms.get_qe_scf_atomic_positions(),
                'kpoints': self.input.wfn.get_kgrid_dft(),
            },
            'kpoints_type': 'crystal',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        input_pw2bgw_dict: dict = {
            'namelists': {
                'input_pw2bgw': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'real_or_complex': '2',
                    'wfng_flag': '.true.',
                    'wfng_file': "'WFN_coo'",
                    'wfng_kgrid': '.true.',
                    'wfng_nk1': self.input_dict['wfn']['kdim'][0],
                    'wfng_nk2': self.input_dict['wfn']['kdim'][1],
                    'wfng_nk3': self.input_dict['wfn']['kdim'][2],
                    'wfng_dk1': 0.0,
                    'wfng_dk2': 0.0,
                    'wfng_dk3': 0.0,
                    'rhog_flag': '.true.',
                    'rhog_file': "'RHO'",
                    'vxcg_flag': '.true.',
                    'vxcg_file': "'VXC'",
                    'vscg_flag': '.true.',
                    'vscg_file': "'VSC'",
                    'vkbg_flag': '.true.',
                    'vkbg_file': "'VKB'",
                }
            }
        }
        input_parabands_dict: dict = {
            'maps': {
                'input_wfn_file': 'WFN_coo',
                'output_wfn_file': 'WFN_parabands.h5',
                'vsc_file': 'VSC',
                'vkb_file': 'VKB',
                'number_bands': self.input_dict['wfn']['num_parabands_cond_bands'] + self.input_dict['total_valence_bands'],

            },
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_wfn_dict['namelists']['system']['noncolin'] = True
            input_wfn_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['wfn']['args']
        args_type = self.input_dict['wfn']['args_type']
        input_wfn_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_wfn_dict
        )
        #pw2bgw.
        args_dict = self.input_dict['wfn']['pw2bgw_args']
        args_type = self.input_dict['wfn']['pw2bgw_args_type']
        input_pw2bgw_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_pw2bgw_dict
        )
        #parabands.
        args_dict = self.input_dict['wfn']['parabands_args']
        args_type = self.input_dict['wfn']['parabands_args_type']
        input_parabands_dict = BgwInputFile.update_dict(
            args_dict,
            args_type,
            input_parabands_dict
        )

        # Get string. 
        self.input_wfn: str = QePwInputFile(input_wfn_dict, self.input_dict).get_input_str()
        self.input_pw2bgw: str = QePwInputFile.write_general(input_pw2bgw_dict)
        self.input_parabands: str = BgwInputFile.write_general(input_parabands_dict)

    def set_jobs_str(self):
        self.job_wfn = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < wfn.in &> wfn.in.out 

cp ./tmp/struct.xml ./wfn.xml
'''

        self.job_wfn_pw2bgw = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_pw2bgw_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_pw2bgw_info)}pw2bgw.x -pd .true. < wfn_pw2bgw.in &> wfn_pw2bgw.in.out 
cp ./tmp/WFN_coo ./
cp ./tmp/RHO ./
cp ./tmp/VXC ./
cp ./tmp/VSC ./
cp ./tmp/VKB ./
'''
        
        self.job_parabands = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_parabands_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_parabands_info)}parabands.cplx.x &> parabands.inp.out 
'''

        self.jobs = [
            './job_wfn.sh',
            './job_wfn_pw2bgw.sh',
            './job_parabands.sh',
        ]

    def create(self):
        write_str_2_f(f'wfn.in', self.input_wfn)
        write_str_2_f(f'job_wfn.sh', self.job_wfn)
        write_str_2_f(f'wfn_pw2bgw.in', self.input_pw2bgw)
        write_str_2_f(f'job_wfn_pw2bgw.sh', self.job_wfn_pw2bgw)
        write_str_2_f(f'parabands.inp', self.input_parabands)
        write_str_2_f(f'job_parabands.sh', self.job_parabands)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
        os.system('rm -rf ./wfn.xml')
        os.system('rm -rf ./WFN_coo')
        os.system('rm -rf ./WFN_parabands.h5')
        os.system('rm -rf ./RHO')
        os.system('rm -rf ./VXC')
        os.system('rm -rf ./VSC')
        os.system('rm -rf ./VKB')
        os.system('rm -rf wfn.in.out')
        os.system('rm -rf wfn_pw2bgw.in.out')
        os.system('rm -rf parabands.inp.out')

class WfnqJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.job_pw2bgw_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['wfnq']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['wfnq']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['wfnq']['job_info'])

        # pw2bgw.
        if isinstance(self.input_dict['wfnq']['job_pw2bgw_info'], str):
            self.job_pw2bgw_info = JobProcDesc.from_job_id(
                self.input_dict['wfnq']['job_pw2bgw_info'],
                self.input_dict,
            )
        else:
            self.job_pw2bgw_info = JobProcDesc(**self.input_dict['wfnq']['job_pw2bgw_info'])

    def set_inputs_str(self):
        #Base. 
        input_wfnq_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'bands',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'nbnd': self.input_dict['wfnq']['num_cond_bands'] + self.input_dict['total_valence_bands'],
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                },
                'electrons': {},
                'ions': {},
                'cell': {},
            },
            'blocks': {
                'atomic_species': self.input.atoms.get_qe_scf_atomic_species(),
                'cell_parameters': self.input.atoms.get_qe_scf_cell(),
                'atomic_positions': self.input.atoms.get_qe_scf_atomic_positions(),
                'kpoints': self.input.wfnq.get_kgrid_dft(),
            },
            'kpoints_type': 'crystal',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        input_pw2bgw_dict: dict = {
            'namelists': {
                'input_pw2bgw': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'real_or_complex': '2',
                    'wfng_flag': '.true.',
                    'wfng_file': "'WFNq_coo'",
                    'wfng_kgrid': '.true.',
                    'wfng_nk1': self.input_dict['wfnq']['kdim'][0],
                    'wfng_nk2': self.input_dict['wfnq']['kdim'][1],
                    'wfng_nk3': self.input_dict['wfnq']['kdim'][2],
                    'wfng_dk1': self.input_dict['wfnq']['qshift'][0]*self.input_dict['wfnq']['kdim'][0],
                    'wfng_dk2': self.input_dict['wfnq']['qshift'][1]*self.input_dict['wfnq']['kdim'][1],
                    'wfng_dk3': self.input_dict['wfnq']['qshift'][2]*self.input_dict['wfnq']['kdim'][2],
                }
            }
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_wfnq_dict['namelists']['system']['noncolin'] = True
            input_wfnq_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['wfnq']['args']
        args_type = self.input_dict['wfnq']['args_type']
        input_wfnq_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_wfnq_dict
        )
        #pw2bgw.
        args_dict = self.input_dict['wfnq']['pw2bgw_args']
        args_type = self.input_dict['wfnq']['pw2bgw_args_type']
        input_pw2bgw_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_pw2bgw_dict
        )

        # Get string. 
        self.input_wfnq: str = QePwInputFile(input_wfnq_dict, self.input_dict).get_input_str()
        self.input_pw2bgw: str = QePwInputFile.write_general(input_pw2bgw_dict)

    def set_jobs_str(self):
        self.job_wfnq = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < wfnq.in &> wfnq.in.out 
'''

        self.job_wfnq_pw2bgw = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_pw2bgw_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_pw2bgw_info)}pw2bgw.x -pd .true. < wfnq_pw2bgw.in &> wfnq_pw2bgw.in.out 
cp ./tmp/WFNq_coo ./
cp ./tmp/struct.xml ./wfnq.xml
wfn2hdf.x BIN WFNq_coo WFNq_coo.h5 
'''

        self.jobs = [
            './job_wfnq.sh',
            './job_wfnq_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f(f'wfnq.in', self.input_wfnq)
        write_str_2_f(f'job_wfnq.sh', self.job_wfnq)
        write_str_2_f(f'wfnq_pw2bgw.in', self.input_pw2bgw)
        write_str_2_f(f'job_wfnq_pw2bgw.sh', self.job_wfnq_pw2bgw)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
        os.system('rm -rf ./wfnq.xml')
        os.system('rm -rf WFNq_coo')
        os.system('rm -rf WFNq_coo.h5')
        os.system('rm -rf wfnq.in.out')
        os.system('rm -rf wfnq_pw2bgw.in.out')

class WfnfiJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.job_pw2bgw_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['wfnfi']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['wfnfi']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['wfnfi']['job_info'])

        # pw2bgw.
        if isinstance(self.input_dict['wfnfi']['job_pw2bgw_info'], str):
            self.job_pw2bgw_info = JobProcDesc.from_job_id(
                self.input_dict['wfnfi']['job_pw2bgw_info'],
                self.input_dict,
            )
        else:
            self.job_pw2bgw_info = JobProcDesc(**self.input_dict['wfnfi']['job_pw2bgw_info'])

    def set_inputs_str(self):
        #Base. 
        input_wfnfi_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'bands',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'nbnd': self.input_dict['wfnfi']['num_cond_bands'] + self.input_dict['total_valence_bands'],
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                },
                'electrons': {},
                'ions': {},
                'cell': {},
            },
            'blocks': {
                'atomic_species': self.input.atoms.get_qe_scf_atomic_species(),
                'cell_parameters': self.input.atoms.get_qe_scf_cell(),
                'atomic_positions': self.input.atoms.get_qe_scf_atomic_positions(),
                'kpoints': self.input.wfnfi.get_kgrid_dft(),
            },
            'kpoints_type': 'crystal',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        input_pw2bgw_dict: dict = {
            'namelists': {
                'input_pw2bgw': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'real_or_complex': '2',
                    'wfng_flag': '.true.',
                    'wfng_file': "'WFN_fii'",
                    'wfng_kgrid': '.true.',
                    'wfng_nk1': self.input_dict['wfnfi']['kdim'][0],
                    'wfng_nk2': self.input_dict['wfnfi']['kdim'][1],
                    'wfng_nk3': self.input_dict['wfnfi']['kdim'][2],
                    'wfng_dk1': 0.0,
                    'wfng_dk2': 0.0,
                    'wfng_dk3': 0.0,
                }
            }
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_wfnfi_dict['namelists']['system']['noncolin'] = True
            input_wfnfi_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['wfnfi']['args']
        args_type = self.input_dict['wfnfi']['args_type']
        input_wfnfi_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_wfnfi_dict
        )
        #pw2bgw.
        args_dict = self.input_dict['wfnfi']['pw2bgw_args']
        args_type = self.input_dict['wfnfi']['pw2bgw_args_type']
        input_pw2bgw_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_pw2bgw_dict
        )

        # Get string. 
        self.input_wfnfi: str = QePwInputFile(input_wfnfi_dict, self.input_dict).get_input_str()
        self.input_pw2bgw: str = QePwInputFile.write_general(input_pw2bgw_dict)

    def set_jobs_str(self):
        self.job_wfnfi = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < wfnfi.in &> wfnfi.in.out 

cp ./tmp/struct.xml ./wfnfi.xml
'''

        self.job_wfnfi_pw2bgw = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_pw2bgw_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_pw2bgw_info)}pw2bgw.x -pd .true. < wfnfi_pw2bgw.in &> wfnfi_pw2bgw.in.out 
cp ./tmp/WFN_fii ./
wfn2hdf.x BIN WFN_fii WFN_fii.h5 
'''

        self.jobs = [
            './job_wfnfi.sh',
            './job_wfnfi_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f(f'wfnfi.in', self.input_wfnfi)
        write_str_2_f(f'job_wfnfi.sh', self.job_wfnfi)
        write_str_2_f(f'wfnfi_pw2bgw.in', self.input_pw2bgw)
        write_str_2_f(f'job_wfnfi_pw2bgw.sh', self.job_wfnfi_pw2bgw)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
        os.system('rm -rf ./wfnfi.xml')
        os.system('rm -rf WFN_fii')
        os.system('rm -rf WFN_fii.h5')
        os.system('rm -rf wfnfi.in.out')
        os.system('rm -rf wfnfi_pw2bgw.in.out')

class WfnqfiJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        self.input_dict: dict = self.input.input_dict
        self.scheduler: Scheduler = Scheduler.from_input_dict(self.input_dict)
        self.job_info: JobProcDesc = None
        self.job_pw2bgw_info: JobProcDesc = None
        self.set_job_info()
        self.set_inputs_str()
        self.set_jobs_str()

    def set_job_info(self):
        if isinstance(self.input_dict['wfnqfi']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['wfnqfi']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['wfnqfi']['job_info'])

        # pw2bgw.
        if isinstance(self.input_dict['wfnqfi']['job_pw2bgw_info'], str):
            self.job_pw2bgw_info = JobProcDesc.from_job_id(
                self.input_dict['wfnqfi']['job_pw2bgw_info'],
                self.input_dict,
            )
        else:
            self.job_pw2bgw_info = JobProcDesc(**self.input_dict['wfnqfi']['job_pw2bgw_info'])

    def set_inputs_str(self):
        #Base. 
        input_wfnqfi_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'bands',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat(),
                    'nbnd': self.input_dict['wfnqfi']['num_cond_bands'] + self.input_dict['total_valence_bands'],
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                },
                'electrons': {},
                'ions': {},
                'cell': {},
            },
            'blocks': {
                'atomic_species': self.input.atoms.get_qe_scf_atomic_species(),
                'cell_parameters': self.input.atoms.get_qe_scf_cell(),
                'atomic_positions': self.input.atoms.get_qe_scf_atomic_positions(),
                'kpoints': self.input.wfnqfi.get_kgrid_dft(),
            },
            'kpoints_type': 'crystal',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        input_pw2bgw_dict: dict = {
            'namelists': {
                'input_pw2bgw': {
                    'outdir': "'./tmp'",
                    'prefix': "'struct'",
                    'real_or_complex': '2',
                    'wfng_flag': '.true.',
                    'wfng_file': "'WFNq_fii'",
                    'wfng_kgrid': '.true.',
                    'wfng_nk1': self.input_dict['wfnqfi']['kdim'][0],
                    'wfng_nk2': self.input_dict['wfnqfi']['kdim'][1],
                    'wfng_nk3': self.input_dict['wfnqfi']['kdim'][2],
                    'wfng_dk1': self.input_dict['wfnq']['qshift'][0]*self.input_dict['wfnqfi']['kdim'][0],
                    'wfng_dk2': self.input_dict['wfnq']['qshift'][1]*self.input_dict['wfnqfi']['kdim'][1],
                    'wfng_dk3': self.input_dict['wfnq']['qshift'][2]*self.input_dict['wfnqfi']['kdim'][2],
                }
            }
        }

        # Additions.
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            input_wfnqfi_dict['namelists']['system']['noncolin'] = True
            input_wfnqfi_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['wfnqfi']['args']
        args_type = self.input_dict['wfnqfi']['args_type']
        input_wfnqfi_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_wfnqfi_dict
        )
        #pw2bgw.
        args_dict = self.input_dict['wfnqfi']['pw2bgw_args']
        args_type = self.input_dict['wfnqfi']['pw2bgw_args_type']
        input_pw2bgw_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=input_pw2bgw_dict
        )

        # Get string. 
        self.input_wfnqfi: str = QePwInputFile(input_wfnqfi_dict, self.input_dict).get_input_str()
        self.input_pw2bgw: str = QePwInputFile.write_general(input_pw2bgw_dict)

    def set_jobs_str(self):
        self.job_wfnqfi = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x {self.scheduler.get_sched_mpi_infix(self.job_info)} < wfnqfi.in &> wfnqfi.in.out 
'''

        self.job_wfnqfi_pw2bgw = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_pw2bgw_info)}

{self.scheduler.get_sched_mpi_prefix(self.job_pw2bgw_info)}pw2bgw.x -pd .true. < wfnqfi_pw2bgw.in &> wfnqfi_pw2bgw.in.out 
cp ./tmp/WFNq_fii ./
cp ./tmp/struct.xml ./wfnqfi.xml
wfn2hdf.x BIN WFNq_fii WFNq_fii.h5 
'''

        self.jobs = [
            './job_wfnqfi.sh',
            './job_wfnqfi_pw2bgw.sh',
        ]

    def create(self):
        write_str_2_f(f'wfnqfi.in', self.input_wfnqfi)
        write_str_2_f(f'job_wfnqfi.sh', self.job_wfnqfi)
        write_str_2_f(f'wfnqfi_pw2bgw.in', self.input_pw2bgw)
        write_str_2_f(f'job_wfnqfi_pw2bgw.sh', self.job_wfnqfi_pw2bgw)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

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
        os.system('rm -rf ./wfnqfi.xml')
        os.system('rm -rf WFNq_fii')
        os.system('rm -rf WFNq_fii.h5')
        os.system('rm -rf wfnqfi.in.out')
        os.system('rm -rf wfnqfi_pw2bgw.in.out')

#endregion
