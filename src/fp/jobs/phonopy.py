#region modules
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
from fp.schedulers.scheduler import JobProcDesc, Scheduler
from fp.inputs.qepw import QePwInputFile, IbravType
import glob 
from ase.dft.kpoints import get_special_points
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class PhonopyJob:
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
        if isinstance(self.input_dict['dfpt']['job_info'], str):
            self.job_info = JobProcDesc.from_job_id(
                self.input_dict['dfpt']['job_info'],
                self.input_dict,
            )
        else:
            self.job_info = JobProcDesc(**self.input_dict['dfpt']['job_info'])

    def set_inputs_str(self):
        # Base.
        qdim = self.input_dict['dfpt']['qdim']
        nqpt = qdim[0] * qdim[1] * qdim[2]
        scf_prefix_dict: dict = {
            'namelists': {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos',
                    'calculation': 'scf',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': IbravType(self.input_dict).get_idx(),
                    'ntyp': self.input.atoms.get_ntyp(),
                    'nat': self.input.atoms.get_nat() * nqpt,
                    'ecutwfc': self.input_dict['scf']['cutoff'],
                },
                'electrons': {},
                'ions': {},
                'cell': {},
            },
            'blocks': {},
            'kpoints_type': 'automatic',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }
        scf_suffix_dict: dict = {
            'namelists': {},
            'blocks': {
                'kpoints': [1, 1, 1] 
            },
            'kpoints_type': 'automatic',   # Options are 'automatic', 'crystal' and 'crystal_b'. 
            'cell_units': self.input_dict['atoms']['write_cell_units'],
            'position_units': self.input_dict['atoms']['write_position_units'],
        }

        # Update. 
        #spinorbit.
        if self.input_dict['scf']['is_spinorbit']:
            scf_prefix_dict['namelists']['system']['noncolin'] = True
            scf_prefix_dict['namelists']['system']['lspinorb'] = True
        #override or extra. 
        args_dict = self.input_dict['scf']['args']
        args_type = self.input_dict['scf']['args_type']
        scf_prefix_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=scf_prefix_dict
        )
        scf_suffix_dict = self.input.update_qe_args_dict(
            args_dict=args_dict,
            args_type=args_type,
            qedict_to_update=scf_suffix_dict
        )
        scf_prefix_dict['blocks'] = {}
        scf_prefix_dict['namelists']['system']['nat'] = self.input.atoms.get_nat() * nqpt
        scf_suffix_dict['namelists'] = {}
        scf_suffix_dict['blocks']['kpoints'] = [1, 1, 1]


        # Get string. 
        self.phonopy_scf_prefix: str = QePwInputFile(scf_prefix_dict, self.input_dict).get_input_str()
        self.phonopy_scf_suffix: str = '\n' + QePwInputFile(scf_suffix_dict, self.input_dict).get_input_str()

    def set_jobs_str(self):
        self.jobs = [
            './job_phonopy.sh',
        ]

    def get_phonopy_bandpath(self) -> str:
        sc_map = get_special_points(self.input.atoms.atoms.cell)
        sc_labels = self.input_dict['path']['special_points']
        npoints_segment = self.input_dict['path']['npoints_segment']
        output = ''

        output += '--band=" '

        for sc_label in sc_labels:
            for col in sc_map[sc_label]:
                output += f' {col:15.10f} '

        output += ' " '

        output += f' --band-points={npoints_segment} '

        return output

    def create_phonopy_files(self):
        # Create supercell files.
        qdim = self.input_dict['dfpt']['qdim']
        os.system(f'phonopy --qe -d --dim="{qdim[0]} {qdim[1]} {qdim[2]}" -c scf.in')
        sc_files = glob.glob('supercell-*')
        for sc_file in sc_files:
            os.system(f'cat phonopy_scf_prefix {sc_file} phonopy_scf_suffix >| phonopy-{sc_file}')


        # Create supercell job.
        phonopy_bandpath_str: str = self.get_phonopy_bandpath()
        files = glob.glob('phonopy-supercell-*')
        start_idx = 0
        stop_idx = len(files)
        debug_str: str = '\n'
        files_bashvar_str: str = '\nfiles=('
        files_args_str: str = ''
        for file_idx, file in enumerate(files): 
            files_bashvar_str += f'"{file}" '
            files_args_str += f' {file}.out '
            debug_str += f'#idx: {file_idx}, filename: {file}\n'
        files_bashvar_str += ')\n\n'
        debug_str += '\n\n'
        file_variable = '${files[$i]}'

        self.job_phonopy = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

{debug_str}

{files_bashvar_str}

start={start_idx}
stop={stop_idx-1}

for (( i=$start; i<=$stop; i++ )); do
{self.scheduler.get_sched_mpi_prefix(self.job_info)}pw.x < {file_variable} &> {file_variable}.out
done

# Post processing. This should create FORCE_SETS and phonon bands. 
phonopy -f {files_args_str}
phonopy --qe -c scf.in {phonopy_bandpath_str} --dim="{qdim[0]} {qdim[1]} {qdim[2]}"
'''

    def create(self):
        write_str_2_f('phonopy_scf_prefix', self.phonopy_scf_prefix)
        write_str_2_f('phonopy_scf_suffix', self.phonopy_scf_suffix)
        self.create_phonopy_files()
        write_str_2_f('job_phonopy.sh', self.job_phonopy)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'job_phonopy.sh',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        os.system('rm -rf phonopy_scf_prefix')
        os.system('rm -rf phonopy_disp.yaml')
        os.system('rm -rf phonopy_scf_suffix')
        os.system('rm -rf phonopy-supercell*')
        os.system('rm -rf supercell*')
        os.system('rm -rf job_phonopy.sh')

        os.system('rm -rf FORCE_SETS')
        os.system('rm -rf phonopy.yaml')
        os.system('rm -rf band.yaml')

#endregions
