#region: Modules.
from ase import Atoms 
import numpy as np 
from ase.data import chemical_symbols, atomic_masses, atomic_numbers
from ase.io import read 
from ase.units import Angstrom, Bohr
import pandas as pd 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
    
class BravaisLattice:
    def __init__(self, ibrav: str, **kwargs):
        self.ibrav: str = ibrav

        for key, value in kwargs.items():
            setattr(self, key, value)

    def get_lattice_vectors(self):
        if self.ibrav=='sc':
            return np.diag([self.A, self.A, self.A])
        elif self.ibrav=='fcc':
            return np.array([
                [-0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
                [-0.5, 0.5, 0.0],
            ])*self.A
        elif self.ibrav=='bcc':
            return np.array([
                [0.5, 0.5, 0.5],
                [-0.5, 0.5, 0.5],
                [-0.5, -0.5, 0.5],
            ])*self.A
        elif self.ibrav=='tetra':
            return np.diag([self.A, self.A, self.C])
        elif self.ibrav=='ortho':
            return np.diag([self.A, self.B, self.C])
        elif self.ibrav=='hex':
            return np.array([
                [self.A, 0.0, 0.0],
                [-0.5*self.A, 0.5*np.sqrt(3)*self.A, 0.0],
                [0.0, 0.0, self.C],
            ])
        else:
            return np.eye(3)

class CellUpdate:
    def __init__(self, input_dict: dict):
        self.input_dict: dict = input_dict

    def get_cell_array(self):
        # Read from input_dict.
        cell_factor = None
        if (cell_list := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('cell', {})
            .get('units'))  is not None:
            cell_units = self.input_dict['atoms']['read']['cell']['units']
            if cell_units=='angstrom':
                cell_factor = Angstrom
            elif cell_units=='bohr':
                cell_factor = Bohr
            elif cell_units=='alat':
                cell_factor = 1.0

        cell_array = None
        if (cell_list := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('cell', {})
            .get('vectors'))  is not None:
            if cell_units=='angstrom':
                cell_array = np.array(cell_list)*cell_factor
            elif cell_units=='bohr':
               cell_array = np.array(cell_list)*cell_factor
            elif cell_units=='alat':
                cell_array = BravaisLattice(
                    **self.input_dict['atoms']['read']['cell']['alat_info']
                ).get_lattice_vectors()*cell_factor

        # Read from file.
        if (cell_file := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('cell', {})
            .get('file'))  is not None:
            if len(open(cell_file).read())!=0: 
                with open(cell_file, 'r') as f: line = f.readline().strip()
                cell_array = np.loadtxt(cell_file, skiprows=1)*cell_factor

        return cell_array

class PositionsUpdate:
    def __init__(self, input_dict: dict):
        self.input_dict: dict = input_dict

    def get_positions_array(self, cell_array: np.ndarray=None):
        # Read from input_dict.
        position_array = None
        if (pos_list := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('positions', {})
            .get('vectors'))  is not None:
            pos_units = self.input_dict['atoms']['read']['positions']['units']
            
            position_array = []
            for row in pos_list: position_array.append(row[1:])
            position_array = np.array(position_array)

            if pos_units=='angstrom':
                position_array *= Angstrom
            elif pos_units=='bohr':
                position_array *= Bohr
            elif pos_units=='crystal':
                position_array = position_array @ cell_array
            elif pos_units=='alat':
                position_array *= position_array * self.input_dict['atoms']['read']['cell']['alat_info']['A']
        
        # Read from file.
        if (pos_file := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('positions', {})
            .get('file'))  is not None:
            pos_factor = None
            if len(open(pos_file).read())!=0: 
                with open(pos_file, 'r') as f: line = f.readline().strip()
                position_array = pd.read_csv(
                    filepath_or_buffer=pos_file,
                    skiprows=1,
                    sep=' ',
                    header=None,
                ).iloc[:, 1:].to_numpy()

                if line=='angstrom':
                    position_array *= Angstrom
                elif line=='bohr':
                    position_array *= Bohr
                elif line=='crystal':
                    position_array = position_array @ cell_array
                elif line=='alat':
                    position_array *= position_array * self.input_dict['atoms']['read']['cell']['alat_info']['A']

        return position_array

    def get_atomic_numbers(self):
        # Read from input_dict.
        atm_numbers = None
        if (pos_list := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('positions', {})
            .get('vectors'))  is not None:
            atm_numbers = []
            for row in pos_list:
                if isinstance(row[0], str):
                    atm_numbers.append(atomic_numbers[row[0]])
                else:
                    atm_numbers.append(row[0])

        # Read from file.
        if (pos_file := self.input_dict
            .get('atoms', {})
            .get('read', {})
            .get('positions', {})
            .get('file'))  is not None:
            if len(open(pos_file).read())!=0: 
                with open(pos_file, 'r') as f: line = f.readline().strip()
                symbol_list = pd.read_csv(
                    filepath_or_buffer=pos_file,
                    skiprows=1,
                    sep=' ',
                    header=None,
                ).iloc[:, 0].to_list()
                atm_numbers = [ atomic_numbers[item] for item in symbol_list]

        return atm_numbers

class AtomsInput:
    def __init__(
        self,
        input_dict: dict,
    ):
        self.input_dict: dict = input_dict

        # Create atoms.
        self.atoms: Atoms = None
        self.update_atoms()

    def update_atoms(self):
        # Read from atoms file.
        if (atoms_file := self.input_dict.get('atoms', {}).get('file'))  is not None:
            self.atoms = read(atoms_file)

        # Cell.
        cell_array: np.ndarray = None
        cell_array = CellUpdate(self.input_dict).get_cell_array()

        # Positions.
        positions_array: np.ndarray = None
        positions_obj = PositionsUpdate(self.input_dict)
        positions_array = positions_obj.get_positions_array(cell_array=cell_array)

        # Atomic numbers.
        atm_numbers: np.ndarray = None
        atm_numbers = positions_obj.get_atomic_numbers()

        # Create atoms.
        if self.atoms is None:
            self.atoms = Atoms(
                numbers=atm_numbers,
                cell=cell_array,
                positions=positions_array,
                pbc=[True, True, True],
            )
        else:
            if atm_numbers is not None: self.atoms.numbers = atm_numbers
            if positions_array is not None: self.atoms.positions = positions_array
            if cell_array is not None: self.atoms.cell.array = cell_array
        
    def get_ntyp(self):
        return len(np.unique(self.atoms.get_atomic_numbers()))

    def get_nat(self):
        return len(self.atoms)
    
    def get_qe_scf_cell(self, fmt=None):
        output = []

        cell_array = self.atoms.get_cell()

        if self.input_dict['atoms']['write_cell_units']=='bohr': cell_array /= Bohr
        if self.input_dict['atoms']['write_cell_units']=='alat': cell_array /= self.input_dict['atoms']['read']['cell']['alat_info']['A']

        for row in cell_array:
            output.append([
                float(row[0]),
                float(row[1]),
                float(row[2]),
            ])

        # If format is string, convert and return early.
        if fmt=='str':
            output_str = ''
            for row in output:
                for col in row:
                    output_str += f'{col} '
                output_str += '\n'
            return output_str


        return output
    
    def get_qe_scf_atomic_species(self):
        output = []
        
        for atm_num in np.unique(self.atoms.get_atomic_numbers()):
            output.append([
                chemical_symbols[atm_num],
                float(atomic_masses[atm_num]),
                f'{chemical_symbols[atm_num]}.upf'
            ])

        return output 
    
    def get_qe_scf_atomic_positions(self, first_column='symbol', fmt=None):
        output = []

        pos_array = self.atoms.get_scaled_positions() if self.input_dict['atoms']['write_position_units']=='crystal' else self.atoms.get_positions()

        if self.input_dict['atoms']['write_position_units']=='bohr': pos_array /= Bohr
        if self.input_dict['atoms']['write_position_units']=='alat': pos_array /= self.input_dict['atoms']['read']['cell']['alat_info']['A']

        if first_column=='symbol':
            for atm_num, row in zip(self.atoms.get_atomic_numbers(), pos_array):
                output.append([
                    chemical_symbols[atm_num],
                    float(row[0]),
                    float(row[1]),
                    float(row[2])
                ])
        elif first_column=='atom_index':
            _, atom_index = np.unique(self.atoms.get_atomic_numbers(), return_inverse=True)
            atom_index += 1     # 1 based index.
            for atm_num, row in zip(atom_index, pos_array):
                output.append([
                    atm_num,
                    float(row[0]),
                    float(row[1]),
                    float(row[2])
                ])
        else:
            for row in pos_array:
                output.append([
                    float(row[0]),
                    float(row[1]),
                    float(row[2]),
                ])

        if fmt=='str':
            output_str = ''
            for row in output:
                for col in row:
                    output_str += f'{col} '
                output_str += '\n'
            return output_str


        return output

    def get_abacus_atomic_species(self):
        output = []
        
        for atm_num in np.unique(self.atoms.get_atomic_numbers()):
            output.append([
                chemical_symbols[atm_num],
                float(atomic_masses[atm_num]),
                f'{chemical_symbols[atm_num]}.upf'
            ])

        return output 

    def get_abacus_orbitals(self):
        output = []
        
        for atm_num in np.unique(self.atoms.get_atomic_numbers()):
            output.append([
                f'{chemical_symbols[atm_num]}.orb'
            ])

        return output 

    def get_abacus_cell(self):
        # Assuming cell units is Angstrom. Not changing this for abacus.
        output = []

        cell_array = self.atoms.get_cell()

        for row in cell_array:
            output.append([
                float(row[0]),
                float(row[1]),
                float(row[2]),
            ])

        return output

    def get_abacus_atomic_positions(self):
        output = []

        pos_array = self.atoms.get_positions()
        atm_numbers = self.atoms.get_atomic_numbers()

        for atm_number, coord in zip(atm_numbers, pos_array):
            output.append([
                chemical_symbols[atm_number],
                float(coord[0]),
                float(coord[1]),
                float(coord[2]),
            ])

        return output

    def get_wan_cell(self):
        output = []
        output.append(['Ang'])
        
        for row in self.atoms.cell:
            output.append(row)

        return output 

    def get_wan_atomic_positions(self):
        output = []
        output.append(['Ang'])

        for atom, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
            output.append([atom, *pos])

        return output 

#endregion
