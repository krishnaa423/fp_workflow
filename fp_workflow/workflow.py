#region: Modules.
from ase import Atoms
from ase.io import read, write 
from ase.build import make_supercell
from ase.data import atomic_masses, atomic_numbers, chemical_symbols
import numpy as np 
import os 
import argparse
import time 
import subprocess
import glob 
#endregion

#region: Variables.

# TODO: Single wfn class and k-points as inputs. Additional steps in the workflow (phbands, xctph, esf). Test on Si. 
# TODO: Convergence. Test on Si. And on SiO2 too. 

total_time:float = 0.0
#endregion

#region: Functions.
def create_atoms():
    
    symbols = ['Si', 'Si']
    A = 5.43
    atoms = Atoms(
        numbers=[atomic_numbers[symbol] for symbol in symbols],
        cell=np.array([
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ])*A,
        scaled_positions=np.array([
            [0.00, 0.00, 0.00],
            [0.25, 0.25, 0.25],
        ]),
        pbc=[True, True, True],
    )
    # symbols = [
    #     'Cs',
    #     'Cs',
    #     'Ag',
    #     'In',
    #     'Cl',
    #     'Cl',
    #     'Cl',
    #     'Cl',
    #     'Cl',
    #     'Cl',
    # ]
    # A = 10.48
    # atoms = Atoms(
    #     numbers=[atomic_numbers[symbol] for symbol in symbols],
    #     cell=np.array([
    #         [0.5, 0.5, 0.0],
    #         [0.0, 0.5, 0.5],
    #         [0.5, 0.0, 0.5],
    #     ])*A,
    #     scaled_positions=np.array([
    #         [0.25, 0.25, 0.25],
    #         [0.75, 0.75, 0.75],
    #         [0.50, 0.50, 0.50],
    #         [0.00, 0.00, 0.00],
    #         [0.25, 0.25, 0.75],
    #         [0.25, 0.75, 0.25],
    #         [0.75, 0.25, 0.25],
    #         [0.75, 0.75, 0.25],
    #         [0.75, 0.25, 0.75],
    #         [0.25, 0.75, 0.75],
    #     ]),
    #     pbc=[True, True, True],
    # )
    # write('struct_ase_uc.cif', [atoms])        # Write the unit cell cif file.

    # # Make supercell. 
    # atoms = make_supercell(atoms, P=np.diag(np.array([2, 2, 2])))
    # write('struct_ase_sc.cif', [atoms])        # Write the supercell cell cif file.
    # atoms.rattle()
    # write('struct_ase_sc_rattled.cif', [atoms])        # Write the supercell cell cif file.
    
    
    
    return atoms 

def create_input_dict():
    
    # scheduler = Summit()
    scheduler = WSL()
    
    atoms = AtomsInput(
        atoms = create_atoms(),
    )
    
    scf = ScfInput(
        kgrid=(2, 2, 2),
        ecutwfc=20.0,
        job_proc_desc=JobProcDesc()
    )
    
    dfpt = DfptInput(
        qgrid= (2, 2, 2),
        conv_threshold='1.0d-10',
        job_proc_desc=JobProcDesc(),
    )
    
    wfn = WfnGeneralInput(
        atoms=atoms,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.000),
        is_sym=True,
        bands=101,
        job_proc_desc=JobProcDesc(),
        parabands_bands=201,
    )
    
    wfnq = WfnGeneralInput(
        atoms=atoms,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.001),
        is_sym=True,
        bands=30,
        job_proc_desc=JobProcDesc(),
    )
    
    wfnfi = WfnGeneralInput(
        atoms=atoms,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.0),
        is_sym=False,
        bands=30,
        job_proc_desc=JobProcDesc(),
    )
    
    epw = EpwInput(
        kgrid_coarse=(2, 2, 2),
        qgrid_coarse=(2, 2, 2),
        kgrid_fine=(2, 2, 2),
        qgrid_fine=(2, 2, 2),
        bands=30,
        job_proc_desc=JobProcDesc(),
    )
    
    wfnqfi = WfnGeneralInput(
        atoms=atoms,
        kdim=(2, 2, 2),
        qshift=(0.0, 0.0, 0.001),
        is_sym=False,
        bands=30,
        job_proc_desc=JobProcDesc(),
    )
    
    epsilon = EpsilonInput(
        wfn_input=wfn,
        bands=200,
        cutoff=10.0,
        job_proc_desc=JobProcDesc()
    )
    
    sigma = SigmaInput(
        wfn_input=wfn,
        bands=200,
        band_min=1,
        band_max=30,
        cutoff=10.0,
        job_proc_desc=JobProcDesc(),
    )
    
    kernel = KernelInput(
        val_bands_coarse=4,
        cond_bands_coarse=6,
        job_proc_desc=JobProcDesc(),
    )
    
    absorption = AbsInput(
        val_bands_coarse=4,
        cond_bands_coarse=6,
        val_bands_fine=4,
        cond_bands_fine=6,
        job_proc_desc=JobProcDesc(),
    )
    
    plotxct = PlotxctInput(
        hole_position=[1, 1, 1],
        supercell_size=[2, 2, 2],
        state=1,
        job_proc_desc=JobProcDesc(),
    )
    
    input = Input(
        scheduler=scheduler,
        atoms=atoms,
        scf=scf,
        dfpt=dfpt,
        wfn=wfn,
        wfnq=wfnq,
        wfnfi=wfnfi,
        epw=epw,
        wfnqfi=wfnqfi,
        epsilon=epsilon,
        sigma=sigma,
        kernel=kernel,
        absorption=absorption,
        plotxct=plotxct,
    )
    # For some reason I need to include these lines below. Don't know why for now. 
    input.scheduler=scheduler
    input.atoms = atoms 
    
    return input 

def write_string_to_file(filename, string):
    with open(filename, 'w') as f:
        f.write(string)
        
def create_workflow(input):
    
    list_of_steps = [
        Scf(input=input),
        Dfpt(input=input),
        Wfn(input=input),
        Wfnq(input=input),
        Wfnfi(input=input),
        Epw(input=input),
        Wfnqfi(input=input),
        Epsilon(input=input),
        Sigma(input=input),
        Kernel(input=input),
        Absorption(input=input),
        PlotXct(input=input),
        # Esf(input=input),
        # Step(input=input),
        # SaveStep(input=input),
        # Update(input=input),
    ]
    
    for step in list_of_steps:
        step.create_files()

    # Add exec permission to scripts. 
    for file in glob.glob('*.sh'):
        os.system(f'chmod u+x {file}') 

def run_and_wait_command(script, input):
    '''
    Run each script and write out some logging info. 
    '''
    global total_time

    start_time = time.time()
    print(f'Starting {script}.', flush=True)
    ps_result = subprocess.run(f'{input.scheduler.job_prefix}{script}')
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    total_time += elapsed_time

    if ps_result.returncode == 0:  # Success.
        print(f'Done with {script} in {elapsed_time:15.10f} seconds.\n\n', flush=True)
    else:               # Fail.
        print(f'Error finishing: {script}. Exited with code {ps_result.returncode}. Time elapsed is {elapsed_time:15.10f} seconds.\n\n', flush=True)
        print(f'Total time for workflow run in {total_time:15.10f} seconds.\n', flush=True)
        os._exit(ps_result.returncode)
        
def run_workflow(input):
    
    list_of_steps = [
        Scf(input=input),
        Dfpt(input=input),
        Wfn(input=input),
        Wfnq(input=input),
        Wfnfi(input=input),
        Epw(input=input),
        Wfnqfi(input=input),
        Epsilon(input=input),
        Sigma(input=input),
        Kernel(input=input),
        Absorption(input=input),
        PlotXct(input=input),
        # Esf(input=input),
        # Step(input=input),
        # SaveStep(input=input),
        # Update(input=input),
    ]
    
    for step in list_of_steps:
        step.run()
        
    # Write the total workflow run time. 
    print(f'Done whole worflow in {total_time:15.10f} seconds.\n\n', flush=True)

def remove_workflow(input):
    
    list_of_steps = [
        Scf(input=input),
        Dfpt(input=input),
        Wfn(input=input),
        Wfnq(input=input),
        Wfnfi(input=input),
        Epw(input=input),
        Wfnqfi(input=input),
        Epsilon(input=input),
        Sigma(input=input),
        Kernel(input=input),
        Absorption(input=input),
        PlotXct(input=input),
        # Esf(input=input),
        # Step(input=input),
        # SaveStep(input=input),
        # Update(input=input),
    ]
    
    for step in list_of_steps:
        step.delete_files()
        
    os.system('rm -rf *.cif')

def save_step(folder_name):
    list_of_steps = [
        Scf(input=input),
        Dfpt(input=input),
        Wfn(input=input),
        Wfnq(input=input),
        Wfnfi(input=input),
        Epw(input=input),
        Wfnqfi(input=input),
        Epsilon(input=input),
        Sigma(input=input),
        Kernel(input=input),
        Absorption(input=input),
        PlotXct(input=input),
        Esf(input=input),
        Step(input=input),
    ]
    
    for step in list_of_steps:
        step.save(folder_name)

def main():
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Script to manage workflows")

    # Add arguments
    parser.add_argument('--create', action='store_true', help="create workflow")
    parser.add_argument('--run', action='store_true', help="run workflow")
    parser.add_argument('--remove', action='store_true', help="remove workflow")

    # Parse the arguments
    args = parser.parse_args()
    
    # Perform actions based on arguments. 
    input = create_input_dict()
    if args.create:
        create_workflow(input=input)

    if args.run:
        run_workflow(input=input)

    if args.remove:
        remove_workflow(input=input) 
#endregion

#region: Classes.
class JobProcDesc:
    def __init__(self, time=None, tasks=None, nk=None, ni=None):
        self.time = time 
        self.tasks = tasks 
        self.nk = nk
        self.ni = ni 

class Scheduler:
    def __init__(self):
        self.job_prefix = None
        self.parallelprefix = None 
        self.jobheader = None 
    
    def get_job_header(self, job_proc_desc):
        pass 
    
    def get_parallel_prefix(self, job_proc_desc):
        pass 
    
    def get_job_prefix(self):
        pass 

class Summit(Scheduler):
    def __init__(self):
        self.job_prefix = 'bsub '
        self.parallelprefix = 'jsrun -n240 -a1 -c7 -g1 -bpacked:7 -EOMP_NUM_THREADS=28 '
        self.jobheader = f'''#BSUB -P cph156
#BSUB -q batch
#BSUB -J struct_job
#BSUB -nnodes 40
#BSUB -W 01:00
        '''

class WSL(Scheduler):
    def __init__(self):
        self.job_prefix = ''
        self.parallelprefix = ''
        self.jobheader = f''''''

class AtomsInput:
    def __init__(
        self,
        atoms,
    ):
        self.atoms:Atoms = atoms  
        
    def get_ntyp(self):
        return len(np.unique(self.atoms.get_atomic_numbers()))
    
    def get_nat(self):
        return self.atoms.get_number_of_atoms()

    def get_scf_cell(self):
        output = ''
        for row in self.atoms.get_cell():
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        return output 
    
    def get_scf_atomic_species(self):
        output = ''
        
        for atm_num in np.unique(self.atoms.get_atomic_numbers()):
            output += f'{chemical_symbols[atm_num]} {atomic_masses[atm_num]} {chemical_symbols[atm_num]}_ONCV_PBE_sr.upf\n'
        return output 

    def get_scf_atomic_positions(self, first_column='symbol'):
        output = ''
        
        if first_column=='symbol':
            for atm_num, row in zip(self.atoms.get_atomic_numbers(), self.atoms.get_positions()):
                output += f'{chemical_symbols[atm_num]} {row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        
        if first_column=='atom_index':
            _, atom_index = np.unique(self.atoms.get_atomic_numbers(), return_inverse=True)
            atom_index += 1     # 1 based index.
            for atm_num, row in zip(atom_index, self.atoms.get_positions()):
                output += f'{atm_num} {row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        return output 

class ScfInput:
    def __init__(
        self,
        kgrid,
        ecutwfc, 
        job_proc_desc,
    ):
        self.kgrid:np.ndarray = np.array(kgrid) 
        self.ecutwfc:float = ecutwfc
        self.job_proc_desc: JobProcDesc = job_proc_desc
        
    def get_kgrid(self):
        output = ''
        output += 'K_POINTS automatic\n'
        output += f'{int(self.kgrid[0])} {int(self.kgrid[1])} {int(self.kgrid[2])} 0 0 0\n'
        
        return output 

class DfptInput:
    def __init__(
        self,
        qgrid,
        conv_threshold,
        job_proc_desc,
    ):
        self.qgrid: np.array = np.array(qgrid)
        self.conv_threshold: float = conv_threshold
        self.job_proc_desc: JobProcDesc = job_proc_desc

class PhbandsInput:
    pass 

class DosInput:
    pass 

class PdosInput:
    pass 

class KpdosInput:
    pass 

class WannierInput:
    pass 

class WfnGeneralInput:
    def __init__(
        self,
        atoms,
        kdim,
        qshift,
        is_sym,
        bands, 
        job_proc_desc,
        parabands_bands=None,
    ):      
        self.atoms:AtomsInput = atoms 
        self.kdim = kdim
        self.qshift = qshift 
        self.is_sym:bool = is_sym
        self.bands = bands  
        self.job_proc_desc: JobProcDesc = job_proc_desc
        self.parabands_bands = parabands_bands 
        
        self.kgrid:np.ndarray = None 
        self.create_kgrid() 
    
    
    def create_kgrid(self):
        if self.is_sym==True: # Use kgrid.x
            with open('kgrid.inp', 'w') as f:
                f.write(f'{self.kdim[0]:15.10f} {self.kdim[1]:15.10f} {self.kdim[2]:15.10f}\n')     
                f.write(f'0.0 0.0 0.0\n')     
                f.write(f'{self.qshift[0]:15.10f} {self.qshift[1]:15.10f} {self.qshift[2]:15.10f}\n')
                f.write(f'{self.atoms.get_scf_cell()}')
                f.write(f'{self.atoms.get_nat()}\n')
                f.write(f'{self.atoms.get_scf_atomic_positions(first_column="atom_index")}')
                f.write(f'0 0 0\n')
                f.write(f'.false.\n')
                f.write(f'.false.\n')
                f.write(f'.false.\n')
            
            os.system('kgrid.x kgrid.inp kgrid.log kgrid.out')
            
            self.kgrid = np.loadtxt('kgrid.log', skiprows=2)
            
            if self.kgrid.ndim == 1 : self.kgrid = self.kgrid.reshape(1, self.kgrid.size)
        else: # Use kmesh.pl 
            os.system(f'kmesh.pl {self.kdim[0]:15.10f} {self.kdim[1]:15.10f} {self.kdim[2]:15.10f} > kgrid.log')
            
            self.kgrid = np.loadtxt('kgrid.log', skiprows=2)
            if self.kgrid.ndim == 1 : self.kgrid = self.kgrid.reshape(1, self.kgrid.size)
            
            # Add the qshift to the non-symmetric grid. 
            self.kgrid[:, 0] += self.qshift[0]
            self.kgrid[:, 1] += self.qshift[1]
            self.kgrid[:, 2] += self.qshift[2]
                    
    def get_kgrid_dft_string(self):
        output = ''
        output += 'K_POINTS crystal\n'
        
        num_kpts = self.kgrid.shape[0]
        output += f'{num_kpts}\n'
        
        for row in self.kgrid:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} {row[3]:15.10f}\n'
        
        return output 
   
    def get_kgrid_eps_string(self, qshift):
        output = ''
        output += 'begin qpoints\n'
        
        # qshift.
        output += f'{qshift[0]:15.10f} {qshift[1]:15.10f} {qshift[2]:15.10f} 1.0 1\n' 
        
        # qgrid. 
        for row_idx, row in enumerate(self.kgrid):
            if row_idx==0: continue 
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} 1.0 0\n'
        
        output += 'end\n'
        
        return output 
    
    def get_kgrid_sig_string(self):
        output = ''
        output += 'begin kpoints\n'
        
        # kgrid. 
        for row in self.kgrid:
            output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} 1.0\n'
        
        output += 'end\n'
        
        return output 
   
class EpwInput:
    def __init__(
        self,
        kgrid_coarse,
        qgrid_coarse,
        kgrid_fine,
        qgrid_fine,
        bands,
        job_proc_desc,
    ): 
        self.kgrid_coarse:np.ndarray = kgrid_coarse
        self.qgrid_coarse:np.ndarray = qgrid_coarse
        self.kgrid_fine:np.ndarray = kgrid_fine
        self.qgrid_fine:np.ndarray = qgrid_fine
        self.bands = bands  
        self.job_proc_desc: JobProcDesc = job_proc_desc
        
class EpsilonInput:
    def __init__(
        self,
        wfn_input,
        bands,
        cutoff,
        job_proc_desc,
    ):
        self.wfn_input: WfnGeneralInput = wfn_input
        self.bands = bands  
        self.cutoff = cutoff
        self.job_proc_desc: JobProcDesc = job_proc_desc
        
    def get_qgrid_str(self, qshift):
        return self.wfn_input.get_kgrid_eps_string(qshift=qshift)
        
class SigmaInput:
    def __init__(
        self,
        wfn_input,
        bands,
        band_min,
        band_max,
        cutoff,
        job_proc_desc,
    ):
        self.wfn_input: WfnGeneralInput = wfn_input
        self.bands = bands 
        self.band_min = band_min  
        self.band_max = band_max 
        self.cutoff = cutoff
        self.job_proc_desc: JobProcDesc = job_proc_desc 
        
    def get_kgrid_str(self):
        return self.wfn_input.get_kgrid_sig_string()

class KernelInput:
    def __init__(
        self,
        val_bands_coarse,
        cond_bands_coarse,
        job_proc_desc,
    ):
        self.val_bands_coarse = val_bands_coarse  
        self.cond_bands_coarse = cond_bands_coarse 
        self.job_proc_desc: JobProcDesc = job_proc_desc

class AbsInput:
    def __init__(
        self,
        val_bands_coarse,
        cond_bands_coarse,
        val_bands_fine,
        cond_bands_fine,
        job_proc_desc,
    ):
        self.val_bands_coarse = val_bands_coarse  
        self.cond_bands_coarse = cond_bands_coarse 
        self.val_bands_fine = val_bands_fine 
        self.cond_bands_fine = cond_bands_fine 
        self.job_proc_desc: JobProcDesc = job_proc_desc

class PlotxctInput:
    def __init__(
        self,
        hole_position,
        supercell_size,
        state,
        job_proc_desc,
    ):
        self.hole_position = hole_position 
        self.supercell_size = supercell_size
        self.state = state 
        self.job_proc_desc: JobProcDesc = job_proc_desc
        
    def get_hole_position_str(self):
        output = f'{self.hole_position[0]:15.10f} {self.hole_position[1]:15.10f} {self.hole_position[2]:15.10f}'
        
        return output 
    
    def get_supercell_size_str(self):
        output = f'{int(self.supercell_size[0])} {int(self.supercell_size[1])} {int(self.supercell_size[2])}'
        
        return output 

class Input: 
    def __init__(
        self,
        scheduler,
        atoms: AtomsInput,
        scf: ScfInput,
        dfpt: DfptInput,
        wfn: WfnGeneralInput,
        wfnq: WfnGeneralInput,
        wfnfi: WfnGeneralInput,
        epw: EpwInput,
        wfnqfi: WfnGeneralInput,
        epsilon: EpsilonInput,
        sigma: SigmaInput,
        kernel: KernelInput,
        absorption: AbsInput,
        plotxct: PlotxctInput, 
    ):
        self.scheduler:Scheduler = scheduler,
        self.atoms:AtomsInput = atoms,
        
        self.scf: ScfInput = scf 
        self.dfpt: DfptInput = dfpt  
        
        self.wfn:WfnGeneralInput = wfn
        self.wfnq:WfnGeneralInput = wfnq
        self.wfnfi:WfnGeneralInput = wfnfi
        self.epw:EpwInput = epw
        self.wfnqfi:WfnGeneralInput = wfnqfi
        
        self.epsilon:EpsilonInput = epsilon
        self.sigma:SigmaInput = sigma
        self.kernel:KernelInput = kernel 
        self.absorption:AbsInput = absorption 
        self.plotxct:PlotxctInput = plotxct 

class Scf:
    def __init__(self, input: Input):
        self.input: Input = input
    
    def create_files(self): 
        
        write_string_to_file(
            'scf.in',
            f'''&CONTROL
outdir='./tmp'
prefix='struct'
pseudo_dir='./ONCVPSP/sg15'
calculation='scf'
tprnfor=.true. 
/

&SYSTEM
ibrav=0
ntyp={self.input.atoms.get_ntyp()}
nat={self.input.atoms.get_nat()}
!nbnd=10
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

{self.input.scf.get_kgrid()}
            '''
        ) 
        
        write_string_to_file(
            'job_scf.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw.x < scf.in &> scf.in.out  
            ''',
        )

    def delete_files(self):
        os.system('rm -rf scf.in')
        os.system('rm -rf job_scf.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf scf.in.out')
        os.system('rm -rf kgrid.inp kgrid.log kgrid.out')

    def run(self):
        run_and_wait_command('./job_scf.sh', self.input)

    def save(self, folder_name): 
        pass 

class Dfpt:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        write_string_to_file(
            'dfpt.in',
            f'''&INPUTPH
outdir='./tmp'
prefix='struct'
ldisp=.true.
nq1={self.input.dfpt.qgrid[0]}
nq2={self.input.dfpt.qgrid[1]}
nq3={self.input.dfpt.qgrid[2]}
tr2_ph={self.input.dfpt.conv_threshold}
fildyn='struct.dyn'
fildvscf='dvscf'
/
            '''
        ) 
        
        write_string_to_file(
            'job_dfpt.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}ph.x < dfpt.in &> dfpt.in.out  
python3 create_save.py
            ''',
        )

        os.system('cp ./utils/create_save.py ./')

    def delete_files(self):
        os.system('rm -rf dfpt.in')
        os.system('rm -rf create_save.py')
        os.system('rm -rf job_dfpt.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf ./save')
        os.system('rm -rf struct.dyn*')
        os.system('rm -rf dfpt.in.out')

    def run(self):
        run_and_wait_command('./job_dfpt.sh', self.input)

    def save(self, folder_name): 
        pass 

class PhBands:
    pass 

class Dos:
    pass 

class Pdos:
    pass 

class Kpdos:
    pass 

class Wannier:
    pass 

class Wfn:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        # wfn. 
        write_string_to_file(
            'wfn.in',
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
        ) 
        
        write_string_to_file(
            'job_wfn.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw.x < wfn.in &> wfn.in.out 
            ''',
        )
        
        # wfn_pw2bgw. 
        write_string_to_file(
            'wfn_pw2bgw.in',
            f'''
&INPUT_PW2BGW
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
            ''',
        )
        
        write_string_to_file(
            'job_wfn_pw2bgw.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw2bgw.x -pd < wfn_pw2bgw.in &> wfn_pw2bgw.in.out 
cp ./tmp/WFN ./
cp ./tmp/RHO ./
cp ./tmp/VXC ./
cp ./tmp/VSC ./
cp ./tmp/VKB ./
            ''',
        )

        # parabands. 
        write_string_to_file(
            'parabands.inp',
            f'''input_wfn_file WFN
output_wfn_file WFN_parabands.h5 

vsc_file VSC 
vkb_file VKB 

number_bands {self.input.wfn.parabands_bands}
            ''',
        )
        
        write_string_to_file(
            'job_parabands.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}parabands.cplx.x &> parabands.inp.out 
            ''',
        )

    def delete_files(self):
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
        os.system('rm -rf ./vxc.dat')
        os.system('rm -rf ./VXC')
        os.system('rm -rf ./VSC')
        os.system('rm -rf ./VKB')
        os.system('rm -rf wfn.in.out')
        os.system('rm -rf wfn_pw2bgw.in.out')
        os.system('rm -rf parabands.inp.out')

    def run(self):
        run_and_wait_command('./job_wfn.sh', self.input)
        run_and_wait_command('./job_wfn_pw2bgw.sh', self.input)
        run_and_wait_command('./job_parabands.sh', self.input)

    def save(self, folder_name): 
        pass 

class Wfnq:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        # wfnq. 
        write_string_to_file(
            'wfnq.in',
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
        ) 
        
        write_string_to_file(
            'job_wfnq.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw.x < wfnq.in &> wfnq.in.out 
            ''',
        )
        
        # wfnq_pw2bgw. 
        write_string_to_file(
            'wfnq_pw2bgw.in',
            f'''
&INPUT_PW2BGW
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
            ''',
        )
        
        write_string_to_file(
            'job_wfnq_pw2bgw.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw2bgw.x -pd < wfnq_pw2bgw.in &> wfnq_pw2bgw.in.out 
cp ./tmp/WFNq ./
wfn2hdf.x BIN WFNq WFNq.h5 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf wfnq.in')
        os.system('rm -rf job_wfnq.sh')
        os.system('rm -rf wfnq_pw2bgw.in')
        os.system('rm -rf job_wfnq_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFNq')
        os.system('rm -rf WFNq.h5')
        os.system('rm -rf wfnq.in.out')
        os.system('rm -rf wfnq_pw2bgw.in.out')

    def run(self):
        run_and_wait_command('./job_wfnq.sh', self.input)
        run_and_wait_command('./job_wfnq_pw2bgw.sh', self.input)

    def save(self, folder_name): 
        pass 

class Wfnfi:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        # wfnq. 
        write_string_to_file(
            'wfnfi.in',
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
        ) 
        
        write_string_to_file(
            'job_wfnfi.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw.x < wfnfi.in &> wfnfi.in.out 
            ''',
        )
        
        # wfnq_pw2bgw. 
        write_string_to_file(
            'wfnfi_pw2bgw.in',
            f'''
&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFN_fi'
wfng_kgrid=.true.
wfng_nk1={self.input.wfnq.kdim[0]}
wfng_nk2={self.input.wfnq.kdim[1]}
wfng_nk3={self.input.wfnq.kdim[2]}
wfng_dk1={self.input.wfnq.qshift[0]*self.input.wfnq.kdim[0]}
wfng_dk2={self.input.wfnq.qshift[1]*self.input.wfnq.kdim[1]}
wfng_dk3={self.input.wfnq.qshift[2]*self.input.wfnq.kdim[2]}
/
            ''',
        )
        
        write_string_to_file(
            'job_wfnfi_pw2bgw.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw2bgw.x -pd < wfnfi_pw2bgw.in &> wfnfi_pw2bgw.in.out 
cp ./tmp/WFN_fi ./
wfn2hdf.x BIN WFN_fi WFN_fi.h5 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf wfnfi.in')
        os.system('rm -rf job_wfnfi.sh')
        os.system('rm -rf wfnfi_pw2bgw.in')
        os.system('rm -rf job_wfnfi_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFN_fi')
        os.system('rm -rf WFN_fi.h5')
        os.system('rm -rf wfnfi.in.out')
        os.system('rm -rf wfnfi_pw2bgw.in.out')

    def run(self):
        run_and_wait_command('./job_wfnfi.sh', self.input)
        run_and_wait_command('./job_wfnfi_pw2bgw.sh', self.input)

    def save(self, folder_name): 
        pass 

class Epw:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        write_string_to_file(
            'epw.in',
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
!bands_skipped='exclude_bands=51:101'

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
        ) 
        
        write_string_to_file(
            'job_epw.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}$START_DIR/Fortran/QuantumEspresso/q-e-cpu/bin/epw.x < epw.in &> epw.in.out 
cp ./tmp/struct_elph* ./
            ''',
        )

    def delete_files(self):
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

    def run(self):
        run_and_wait_command('./job_epw.sh', self.input)

    def save(self, folder_name): 
        pass 

class Wfnqfi:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        # wfnq. 
        write_string_to_file(
            'wfnqfi.in',
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
        ) 
        
        write_string_to_file(
            'job_wfnqfi.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw.x < wfnqfi.in &> wfnqfi.in.out 
            ''',
        )
        
        # wfnq_pw2bgw. 
        write_string_to_file(
            'wfnqfi_pw2bgw.in',
            f'''
&INPUT_PW2BGW
outdir='./tmp'
prefix='struct'
real_or_complex=2
wfng_flag=.true.
wfng_file='WFNq_fi'
wfng_kgrid=.true.
wfng_nk1={self.input.wfnqfi.kdim[0]}
wfng_nk2={self.input.wfnqfi.kdim[1]}
wfng_nk3={self.input.wfnqfi.kdim[2]}
wfng_dk1={self.input.wfnqfi.qshift[0]*self.input.wfnqfi.kdim[0]}
wfng_dk2={self.input.wfnqfi.qshift[1]*self.input.wfnqfi.kdim[1]}
wfng_dk3={self.input.wfnqfi.qshift[2]*self.input.wfnqfi.kdim[2]}
/
            ''',
        )
        
        write_string_to_file(
            'job_wfnqfi_pw2bgw.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}pw2bgw.x -pd < wfnqfi_pw2bgw.in &> wfnqfi_pw2bgw.in.out 
cp ./tmp/WFNq_fi ./
wfn2hdf.x BIN WFNq_fi WFNq_fi.h5 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf wfnqfi.in')
        os.system('rm -rf job_wfnqfi.sh')
        os.system('rm -rf wfnqfi_pw2bgw.in')
        os.system('rm -rf job_wfnqfi_pw2bgw.sh')
        
        os.system('rm -rf ./tmp')
        os.system('rm -rf WFNq_fi')
        os.system('rm -rf WFNq_fi.h5')
        os.system('rm -rf wfnqfi.in.out')
        os.system('rm -rf wfnqfi_pw2bgw.in.out')

    def run(self):
        run_and_wait_command('./job_wfnqfi.sh', self.input)
        run_and_wait_command('./job_wfnqfi_pw2bgw.sh', self.input)

    def save(self, folder_name): 
        pass 

class Epsilon:
    def __init__(self, input):
        self.input:Input = input

    def create_files(self): 
        write_string_to_file(
            'epsilon.inp',
            f'''# Qpoints 
{self.input.epsilon.get_qgrid_str(self.input.wfnq.qshift)}

# Bands
number_bands {self.input.epsilon.bands}
degeneracy_check_override

# G-Cutoff. 
epsilon_cutoff {self.input.epsilon.cutoff}

# Options

# IO. 
use_wfn_hdf5
            '''
        ) 
        
        write_string_to_file(
            'job_epsilon.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

ln -sf ./WFN_parabands.h5 ./WFN.h5 
ln -sf ./WFNq.h5 ./WFNq.h5 
{self.input.scheduler.parallelprefix}epsilon.cplx.x &> epsilon.inp.out 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf epsilon.inp')
        os.system('rm -rf job_epsilon.sh')
        
        os.system('rm -rf ./WFN.h5')
        os.system('rm -rf ./WFNq.h5')
        os.system('rm -rf ./epsmat.h5')
        os.system('rm -rf ./eps0mat.h5')
        os.system('rm -rf epsilon.log')
        os.system('rm -rf chi_converge.dat')
        os.system('rm -rf epsilon.inp.out')

    def run(self):
        run_and_wait_command('./job_epsilon.sh', self.input)

    def save(self, folder_name): 
        pass 

class Sigma:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        write_string_to_file(
            'sigma.inp',
            f'''# kpoints
{self.input.sigma.get_kgrid_str()}
no_symmetries_q_grid
# use_symmetries_q_grid

# Bands.
number_bands {self.input.sigma.bands}
band_index_min {self.input.sigma.band_min}
band_index_max {self.input.sigma.band_max}
degeneracy_check_override

# G-cutoff
screened_coulomb_cutoff {self.input.sigma.cutoff}


# Options
dont_use_vxcdat

# IO
# verbosity 3
use_wfn_hdf5
            '''
        ) 
        
        write_string_to_file(
            'job_sigma.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

ln -sf ./WFN_parabands.h5 ./WFN_inner.h5 
{self.input.scheduler.parallelprefix}sigma.cplx.x &> sigma.inp.out 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf sigma.inp')
        os.system('rm -rf job_sigma.sh')
        
        os.system('rm -rf ./WFN_inner.h5')
        os.system('rm -rf eqp0.dat')
        os.system('rm -rf eqp1.dat')
        os.system('rm -rf sigma_hp.log')
        os.system('rm -rf ch_converge.dat')
        os.system('rm -rf sigma.inp.out')

    def run(self):
        run_and_wait_command('./job_sigma.sh', self.input)

    def save(self, folder_name): 
        pass 

class Sig2Wan:
    pass 

class Kernel:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        write_string_to_file(
            'kernel.inp',
            f'''# Q-points
#exciton_Q_shift 0 Qx Qy Qz
use_symmetries_coarse_grid

# Bands 
number_val_bands {self.input.absorption.val_bands_coarse}
number_cond_bands {self.input.absorption.cond_bands_coarse}
#spinor

# Options
#extended_kernel

# IO. 
use_wfn_hdf5
            '''
        ) 
        
        write_string_to_file(
            'job_kernel.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

ln -sf WFN_parabands.h5 WFN_co.h5
{self.input.scheduler.parallelprefix}kernel.cplx.x &> kernel.inp.out
            ''',
        )

    def delete_files(self):
        os.system('rm -rf kernel.inp')
        os.system('rm -rf job_kernel.sh')
        
        os.system('rm -rf ./WFN_co.h5')
        os.system('rm -rf bsemat.h5')
        os.system('rm -rf kernel.inp.out')

    def run(self):
        run_and_wait_command('./job_kernel.sh', self.input)

    def save(self, folder_name): 
        pass 

class Absorption:
    def __init__(self, input):
        self.input:Input = input

    def create_files(self): 
        write_string_to_file(
            'absorption.inp',
            f'''# Q-points
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
polarization 0.0 0.0 0.001
eqp_co_corrections

# IO
use_wfn_hdf5

# Output
energy_resolution 0.01
write_eigenvectors 10
            '''
        ) 
        
        write_string_to_file(
            'job_absorption.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

ln -sf WFN_parabands.h5 WFN_co.h5 
ln -sf WFN_fi.h5 WFN_fi.h5 
ln -sf WFNq_fi.h5 WFNq_fi.h5 
ln -sf eqp1.dat eqp_co.dat 
{self.input.scheduler.parallelprefix}absorption.cplx.x &> absorption.inp.out 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf absorption.inp')
        os.system('rm -rf job_absorption.sh')
        
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
        os.system('rm -rf bandstructure.dat')
        os.system('rm -rf eigenvectors.h5')
        os.system('rm -rf x.dat')
        os.system('rm -rf epsdiag.dat')
        os.system('rm -rf dtmat')
        os.system('rm -rf vmtxel')
        os.system('rm -rf absorption.inp.out')

    def run(self):
        run_and_wait_command('./job_absorption.sh', self.input)

    def save(self, folder_name): 
        pass 

class PlotXct:
    def __init__(self, input):
        self.input:Input = input
        
    def create_files(self): 
        write_string_to_file(
            'plotxct.inp',
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
        ) 
        
        write_string_to_file(
            'job_plotxct.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct.xsf xsf false abs2 true 
            ''',
        )

    def delete_files(self):
        os.system('rm -rf plotxct.inp')
        os.system('rm -rf job_plotxct.sh')
        
        os.system('rm -rf *.a3Dr')
        os.system('rm -rf plotxct.xsf')
        os.system('rm -rf plotxct.inp.out')

    def run(self):
        run_and_wait_command('./job_plotxct.sh', self.input)

    def save(self, folder_name): 
        pass 

class AbsorptionQ:
    pass 

class XctPh:
    pass 

class Esf:
    def __init__(self, input):
        self.input:Input = input

    def create_files(self): 
            
        write_string_to_file(
            'job_esf.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}python3 esf.py &> esf.out 
            ''',
        )
        
        os.system('cp ./utils/esf.py ./')

    def delete_files(self):
        os.system('rm -rf esf.py')
        os.system('rm -rf job_esf.sh')
        
        os.system('rm -rf esf.h5')
        os.system('rm -rf esf.out')

    def run(self):
        run_and_wait_command('./job_esf.sh', self.input)

    def save(self, folder_name): 
        pass 

class Step:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        write_string_to_file(
            'step.inp',
            f'''
alpha 0.001
max_iter 10

# Max error in force. 
max_error 1.0e-1
            '''
        )
        
        write_string_to_file(
            'job_step.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}python3 step.py &> step.out 
            ''',
        )
        
        os.system('cp ./utils/step.py ./')

    def delete_files(self):
        os.system('rm -rf step.inp')
        os.system('rm -rf step.py')
        os.system('rm -rf job_step.sh')
        
        os.system('rm -rf step.h5')
        os.system('rm -rf step.out')

    def run(self):
        run_and_wait_command('./job_step.sh', self.input)

    def save(self, folder_name): 
        pass 

class SaveStep:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        write_string_to_file(
            'job_savestep.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}python3 savestep.py &> savestep.out 
            ''',
        )
        
        os.system('cp ./utils/savestep.py ./')

    def delete_files(self):
        os.system('rm -rf savestep.py')
        os.system('rm -rf job_savestep.sh')
        
        os.system('rm -rf savestep.out')

    def run(self):
        run_and_wait_command('./job_savestep.sh', self.input)

    def save(self, folder_name): 
        pass 

class Update:
    def __init__(self, input):
        self.input:Input = input
    
    def create_files(self): 
        
        write_string_to_file(
            'job_update.sh',
            f'''#!/bin/bash
{self.input.scheduler.jobheader}

{self.input.scheduler.parallelprefix}python3 update.py &> update.out 
            ''',
        )
        
        os.system('cp ./utils/update.py ./')

    def delete_files(self):
        os.system('rm -rf update.py')
        os.system('rm -rf job_update.sh')
        
        os.system('rm -rf update.out')

    def run(self):
        run_and_wait_command('./job_update.sh', self.input)

    def save(self, folder_name): 
        pass 

#endregion

#region: Main.
if __name__ == "__main__":
    main()
#endregion