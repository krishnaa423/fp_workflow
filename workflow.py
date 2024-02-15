# region: Modules.
import numpy as np 
from ase import Atoms 
from ase.data import chemical_symbols
from ase.build import make_supercell
import json 
from string import Template
# test modules. 
import unittest 
from typing import List, Callable
import os 
import argparse
# endregion

# region: Variables
# endregion

# region: Functions.
def get_str_from_numpy(numpy_array):
    return np.array2string(numpy_array).replace('[', '').replace(']', '')
# endregion

#region: Classes.
class InputJson:
    
    def __init__(self, input_file: str):
        self.input_dict: dict = {
            "struct_alat": '',
            "struct_symbols": '',
            "struct_cell": '',
            "struct_numbers": '',
            "struct_positions": '',
            "struct_supercell": '',
            "struct_kpoints": '',
            "sched_bash_prefix": '',
            "sched_exec": '',
            "sched_launch": '',
            "sched_account": '',
            "sched_constraint": '',
            "sched_qos": '',
            "scf_ecutwfc": '',
            "scf_job_nodes": '',
            "scf_job_tasks": '',
            "scf_job_time": '',
            "wfn_bands": '',
            "wfn_parabands": '',
            "wfn_job_nodes": '',
            "wfn_job_tasks": '',
            "wfn_job_time": '',
            "wfn_job_parabands_nodes": '',
            "wfn_job_parabands_tasks": '',
            "wfn_job_parabands_time": '',
            "wfnq_bands": '',
            "wfnq_job_nodes": '',
            "wfnq_job_tasks": '',
            "wfnq_job_time": '',
            "wfn_fi_bands": '',
            "wfn_fi_ecut": '',
            "wfn_fi_job_nodes": '',
            "wfn_fi_job_tasks": '',
            "wfn_fi_job_time": '',
            "dfpt_cutoff": '',
            "dfpt_job_nodes": '',
            "dfpt_job_tasks": '',
            "dfpt_job_time": '',
            "elph_job_nodes": '',
            "elph_job_tasks": '',
            "elph_job_time": '',
            "elband_kpts": '',
            "elband_job_nodes": '',
            "elband_job_tasks": '',
            "elband_job_time": '',
            "elband_ndiv": '',
            "elband_nkpts_per_div": '',
            "elband_ends": '',
            "fold2bloch_sc": '',
            "fold2bloch_job_nodes": '',
            "fold2bloch_job_tasks": '',
            "fold2bloch_job_time": "01:30:00",
            "wannier_kpts": '',
            "wannier_bands": '',
            "wannier_job_nodes": '',
            "wannier_job_tasks": '',
            "wannier_job_time": '',
            "phband_kpts": '',
            "phband_job_nodes": '',
            "phband_job_tasks": '',
            "phband_job_time": '',
            "phband_job_anaddb_nodes": '',
            "phband_job_anaddb_tasks": '',
            "phband_job_anaddb_time": '',
            "phband_ndiv": '',
            "phband_nkpts_per_div": '',
            "phband_ends": '',
            "eps_bands": '',
            "eps_cutoff": '',
            "eps_job_nodes": '',
            "eps_job_tasks": '',
            "eps_job_time": '',
            "sig_bands": '',
            "sig_job_nodes": '',
            "sig_job_tasks": '',
            "sig_job_time": '',
            "abs_coarse_bands": '',
            "abs_fine_bands": '',
            "abs_job_nodes": '',
            "abs_job_tasks": '',
            "abs_job_time": '',
            "plotxct_state": '',
            "plotxct_hole_position": '',
            "plotxct_job_nodes": '',
            "plotxct_job_tasks": '',
            "plotxct_job_time": '',
            "esf_job_nodes": '',
            "esf_job_tasks": '',
            "esf_job_time": '',
            "step_alpha": '',
            "step_job_nodes": '',
            "step_job_tasks": '',
            "step_job_time": '',
            "update_job_nodes": '',
            "update_job_tasks": '',
            "update_job_time": '',
            "esd_niter": '',
            "esd_iter": ''
        } 

        with open(input_file, 'r') as f: self.input_dict = json.load(f)
        
    def writeJson(self, output_filename: str):
        with open(output_filename, 'w') as f: json.dump(self.input_dict, f)

class Structure:
    
    def __init__(self, **kwargs):
        '''
        Set cell, (positions|scaled_positions), numbers, pbc in the least. 
        '''
        self.cell: np.ndarray = None 
        self.alat: np.ndarray = None 
        self.positions: np.ndarray = None 
        self.scaled_positions: np.ndarray = None 
        self.numbers: np.ndarray = None 
        self.pbc: np.ndarray = None 
        self.ase_atoms: Atoms =  None 
        
        for key, value in kwargs:
            setattr(self, key, value)
        
        if self.positions:
            self.ase_atoms = Atoms(
                cell=self.cell,
                positions=self.positions,
                numbers=self.numbers,
                pbc=self.pbc,
            )
        else:
            self.ase_atoms = Atoms(
                cell=self.cell,
                positions=self.positions,
                numbers=self.numbers,
                pbc=self.pbc,
            )
            
    def getSymbols(self) -> list:
        assert self.ase_atoms
        
        return [chemical_symbols[atomic_number] for atomic_number in self.ase_atoms.get_atomic_numbers().tolist()] 
    
    def getCellStr(self, in_alat_units=False) -> str:
        
        assert self.ase_atoms
        
        cell_array: np.ndarray = self.ase_atoms.get_cell()
        
        # By default it has alat multiplied. 
        if in_alat_units: cell_array /= self.alat 
        
        return np.array2string(cell_array).replace('[', '').replace(']', '')
    
    def getPositionsStr(self, in_alat_units=False) -> str:
        assert self.ase_atoms
        
        cell_array: np.ndarray = self.ase_atoms.get_positions()
        
        # By default it has alat multiplied. 
        if in_alat_units: cell_array /= self.alat 
        
        return np.array2string(cell_array).replace('[', '').replace(']', '')
    
    def getScaledPositionsStr(self) -> str:
        assert self.ase_atoms
        
        cell_array: np.ndarray = self.ase_atoms.get_scaled_positions()
        
        return np.array2string(cell_array).replace('[', '').replace(']', '') 
     
    def genGetSupercell(self, supercell_dims: list | np.ndarray) -> Atoms:
        supercell_dims = np.diag(np.array(supercell_dims))
        
        self.ase_atoms = make_supercell(self.ase_atoms, supercell_dims)
        self.cell = self.ase_atoms.get_cell()
        self.positions = self.ase_atoms.get_positions()
        self.scaled_positions = self.ase_atoms.get_scaled_positions()
        self.numbers = self.ase_atoms.get_atomic_numbers()
        
        return self.ase_atoms

class SchedulerType:
    # Declare some class scoped constants to use 
    INTERACTIVE = 0
    BSUB = 1
    SLURM = 2
                              
class Scheduler:
    
    def __init__(self, input_json: InputJson):
        input_dict: dict = input_json.input_dict
        self.scheduler_exec: str = input_dict['sched_exec'] 
        self.scheduler_prefix: str = input_dict['sched_bash_prefix'] 
        self.account_prefix: str = input_dict['sched_account'] 
        self.queue_prefix: str = input_dict['sched_qos']  
        self.nodes_prefix: str = None 
        self.tasks_prefix: str = None 
        self.time_prefix: str = None 
        self.additional_lines: str = None 
        self.mpi_exec: str = input_dict['sched_exec']  
        self.type = None 
        
        if self.scheduler_exec=='mpirun':
            self.type = SchedulerType.INTERACTIVE
        elif self.scheduler_exec=='jsrun':
            self.type = SchedulerType.BSUB
        elif self.scheduler_exec in ['srun', 'ibrun']:
            self.type = SchedulerType.SLURM
        
    def getJobPrefixExecSuffix(self, queue, nodes, tasks, time) -> str:
        job_prefix: str = '' 
        exec_suffix: str = '' 
        
        if self.type==SchedulerType.INTERACTIVE:
            pass 
        elif self.type==SchedulerType.BSUB:
            pass
        elif self.type==SchedulerType.SLURM:
            pass 
        
        return job_prefix, exec_suffix

class Calculation:
    
    def __init__(self):
        self.job_file_name: str = None  
        self.input_file_name: str = None 
    
    def createFileFromTemplate(self, template: Template, mapping: dict, file_name: str):
        '''
        Generic way to create a mapping from an template and save to a file. 
        '''
        with open(file_name, 'w') as f: f.write(template.substitute(mapping))
    
    def createInputFile(self):
        pass 

    def createJobFile(self):
        pass 
    
    def createFiles(self):
        pass 
    
    def runCalculation(self):
        pass 
    
class Kgrid(Calculation):
    
    def __init__(self, kgrid: tuple | np.ndarray, **kwargs):
        '''
        Override any with kwargs. For example job file name with job_file_name='job_kgrid_wfn.sh'. 
        '''
        super().__init__()
        self.dims: np.ndarray = np.array(kgrid, dtype='i4')        
        self.kshift: np.ndarray = None 
        self.qshift: np.ndarray = None 
        # Set superclass parameters. 
        self.input_file_name = 'kgrid.inp'
        self.job_file_name = 'job_kgrid.sh'
        
        # Override any if necessary.
        for key, value in kwargs: setattr(self, key, value)
        
    def getInputTemplate(self):
        input_file_template = Template(
f'''$kgrid
$kshift
$qshift
$cell_cart_alat
$num_atoms
$pos_cart_alat
0 0 0
.false.
.false.
.false.
'''
        )
        
        return input_file_template
    
    def getJobTemplate(self):
        job_file_template = Template(
f'''$#!/bin/bash
$job_prefix

kgrid.x $exec_suffix $input_file_name kgrid.log kgrid.out
'''
        )
        
        return job_file_template
                
    def createInputFile(self, struct: Structure):
        mapping = {
            'kgrid': get_str_from_numpy(self.dims),
            'kshift': get_str_from_numpy(self.kshift),
            'cell_cart_alat': '',
            'num_atoms': '',
            'pos_cart_alat': '',
        }
        
        self.createFileFromTemplate(
            self.getInputTemplate(),
            mapping,
            self.input_file_name
        )
    
    def createJobFile(self, scheduler: Scheduler, input_json: InputJson):
        # Get job and exec prefix. 
        job_prefix, exec_suffix = scheduler.getJobPrefixExecSuffix(
            input_json.input_dict['sched_qos'],  # queue.
            1,  # nodes.
            1,  # tasks.
            '00:30',  # time.
        )
        
        # Set the mapping. 
        mapping = {
            'job_prefix': job_prefix,
            'exec_prefix': exec_suffix,
        }
        
        # Create the file. 
        self.createFileFromTemplate(
            self.getInputTemplate(),
            mapping,
            self.job_file_name
        )

class Scf(Calculation):
    
    def genQeFiles():
        pass 
    
    def genAbinitFiles():
        pass 

class Nscf(Calculation):
    pass

class Wannier(Calculation):
    pass

class Elbands(Calculation):
    pass

class Fold2Bloch(Calculation):
    pass

class Dfpt(Calculation):
    pass

class Phbands(Calculation):
    pass

class Wfn(Calculation):
    pass

class Wfnq(Calculation):
    pass

class WfnFi(Calculation):
    pass

class WfnqFi(Calculation):
    pass

class Gw(Calculation):
    pass

class Gwbands(Calculation):
    pass

class Bse(Calculation):
    pass

class AbinitAll(Calculation):
    pass 

class CreateWorkflow:
    
    def __init__(self):
        self.calc_list: List[Calculation] = None 
        
    def runWorflow(self):
        for calc in self.calc_list: calc.createFiles()

class RemoveWorkflow:
    
    def __init__(self):
        self.remove_files_list: List[str] = None 
        
    def removeWorkflow(self):
        for filename in self.remove_files_list: os.system(f'rm -rf {filename}')

class RunWorkflow:
    
    def __init__(self):
        self.calc_list: List[Calculation] = None 
        
    def runWorflow(self):
        for calc in self.calc_list: calc.runCalculation()

class Interactive:
    pass 

class CommandLineParser:
    parser = argparse.ArgumentParser(description='Program to manage workflows.')
    
    parser.add_argument('--create')
    parser.add_argument('--remove')
    parser.add_argument('--run')
    
    args = parser.parse_args()
    
    # TODO. 
    if args.create:
        pass 
    elif args.remove:
        pass
    elif args.run:
        pass 


class Bfgs:
    
    def __init__(self, dim: int):
        self.dim: int = dim     # Set the dimension in which the minimization is done.
        self.cur_value: float = None 
        self.cur_gradient: np.ndarray = None 
        self.prev_gradient: np.ndarray = None 
        self.cur_p: np.ndarray = None 
        self.cur_s: np.ndarray = None 
        self.cur_H: np.ndarray = None 
        self.prev_H: np.ndarray = None 
        self.cur_y: np.ndarray = None 
        self.alpha: float = None 
        self.max_iter: int = None 
        self.max_tol: float = None 
        self.converged: bool = False 
        
    def step(self, gradient: np.ndarray):
        self.cur_gradient = gradient
        
        # Obtain direction. 
        
        # Obtain new position. 
        
        # Update y. 
        
        # Update 
    
    def update(self):
        pass
    
    def run(self, function: Callable[[np.ndarray], np.ndarray]):
        pass 


class ReadHdf5:
    
    def __init__(self, file_name):
        self.file_name: str = file_name
        
    def getDatasets(self, ds_names: List[str], slices: List[slice] = None):
        pass 

class ReadNetcdf:
    
    def __init__(self, file_name):
        self.file_name: str = file_name
        
    def getDatasets(self, ds_names: List[str], slices: List[slice] = None):
        pass 

class ReadXml:
    pass 

class WriteHdf5:
    pass 

class DftForce:
    pass 

class Dwfn:
    pass

class Eqp:
    pass 

class Hbse:
    pass

class XctEvec:
    pass 

class Dgw:
    pass 

class P:
    pass

class DK:
    pass 

class Dhbse:
    pass 

class Esf:
    pass 

class XctPh:
    pass 

class Step: 
    pass 

class Update:
    pass 


# region: Test classes
class TestStr2Numpy(unittest.TestCase):
    def testStringConversion(self):
        # Arrange.
        input = np.array([1, 2, 3], dtype='i4') 
        
        # Act. 
        output = get_str_from_numpy(input)
        
        # Assert. 
        self.assertEqual(output, '1 2 3', 'Expected: 1 2 3')
        
class TestBfgs(unittest.TestCase):
    pass 

class TestHbse(unittest.TestCase):
    pass

class TestElph(unittest.TestCase):
    pass 

class TestEsf(unittest.TestCase):
    pass 

class TestXctPh(unittest.TestCase):
    pass 
# endregion

# endregion

# region: Main.
if __name__=='__main__':
    # If testing.
    unittest.main()
    
# endregion