#region modules
from typing import List
from ase.data import atomic_numbers, chemical_symbols
import numpy as np
import os 
from importlib.util import find_spec
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class AbacusInputFile:
    def __init__(
        self,
        abacus_dict: dict,
    ):
        self.abacus_dict: dict = abacus_dict

    @staticmethod
    def from_strings_to_dict(
        input_file: str,
        stru_file: str,
        kpt_file: str,
    ) -> dict:
        pass

    @staticmethod
    def from_dict_to_strings(abacus_dict: dict) -> List[str]:
        pass

    def copy_pseudos_and_orbitals(self):
        # Get directories.
        orb_dir: str = os.path.dirname(find_spec('fp').origin) + '/data/orbitals/abacus'
        pseudo_dir: str = os.path.dirname(find_spec('fp').origin) + '/data/pseudos/abacus'

        # Get unique symbols.
        pos_input_list: list = self.abacus_dict['structure']['atomic_positions']
        atm_symbols: np.ndarray = np.array([ row[0].strip() for row in pos_input_list], dtype='U10')
        unique_atm_symbols = np.unique(atm_symbols)

        # Copy atom pseudos and orbitals.
        os.system('mkdir -p ./abacus_orbitals')
        os.system('mkdir -p ./abacus_pseudos')
        for unique_atm_symbol in unique_atm_symbols:
            os.system(f'cp {orb_dir}/{str(unique_atm_symbol)}.orb ./abacus_orbitals/{str(unique_atm_symbol)}.orb')
            os.system(f'cp {pseudo_dir}/{str(unique_atm_symbol)}.upf ./abacus_pseudos/{str(unique_atm_symbol)}.upf')
        

    def get_input_str(self):
        output = 'INPUT_PARAMETERS\n'

        if self.abacus_dict.get('input') is not None:
            for key, value in self.abacus_dict['input'].items():
                output += f'{key} {value}\n'

        return output

    def get_kpt_str(self):
        #TODO: in case of kpath for dftelbands.
        output = 'K_POINTS\n'

        if self.abacus_dict.get('kpts') is not None:
            kpts = self.abacus_dict['kpts']
            if not isinstance(kpts[0], list):
                output += '0\n'
                output += 'Gamma\n'
                output += f'{kpts[0]} {kpts[1]} {kpts[2]} 0 0 0\n'
            else:   # TODO. kpath for dftelbands.
                pass

        return output

    def get_pseudos(self):
        output = 'ATOMIC_SPECIES\n'

        for row in self.abacus_dict['structure']['atomic_species']:
            for col in row:
                output += f'{col} '
            output += '\n'
        output += '\n'

        return output

    def get_orbitals(self):
        output = 'NUMERICAL_ORBITAL\n'

        for row in self.abacus_dict['structure']['numerical_orbital']:
            for col in row:
                output += f'{col} '
            output += '\n'
        output += '\n'

        return output

    def get_lattice_constant(self):
        output = 'LATTICE_CONSTANT\n'

        for row in self.abacus_dict['structure']['lattice_constant']:
            for col in row:
                output += f'{col} '
            output += '\n'
        output += '\n'

        return output

    def get_cell(self):
        output = 'LATTICE_VECTORS\n'

        for row in self.abacus_dict['structure']['lattice_vectors']:
            for col in row:
                output += f'{col} '
            output += '\n'
        output += '\n'

        return output

    def get_positions(self):
        pos_input_list = self.abacus_dict['structure']['atomic_positions']
        output = 'ATOMIC_POSITIONS\n'
        output += 'Direct\n'

        # Get atomic_numbers list, unique and inverse arrays.
        atm_numbers = np.array(
            [ atomic_numbers[row[0]] for row in pos_input_list],
            dtype='i4',
        )
        pos_array = np.array(
            [ [row[1], row[2], row[3]] for row in pos_input_list],
            dtype='f8',
        )
        unique_atm_numbers = np.unique(
            atm_numbers,
        )

        # Write output.
        for atm_number in unique_atm_numbers:
            output += f'{chemical_symbols[atm_number]}\n'
            output += f'0.0\n' # Magnetic moment I guess.
            
            # Write positions.
            row_idxs = np.where(atm_numbers == atm_number)[0]
            output += f'{len(row_idxs)}\n'
            for row_idx in row_idxs:
                pos = pos_array[row_idx, :]
                output += f'{pos[0]} {pos[1]} {pos[2]} 0 0 0\n'  # Position and then movement. We set the movement to zero. Don't know what it does.

        output += '\n'

        return output

    def get_stru_str(self):
        output = ''

        output += self.get_pseudos()
        output += self.get_orbitals()
        output += self.get_lattice_constant()
        output += self.get_cell()
        output += self.get_positions()

        return output
#endregions
