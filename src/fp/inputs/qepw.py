#region modules
from ase.units import Angstrom, Bohr
import numpy as np
from lark import Lark, Tree
from typing import List
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class NamelistGrammar:
    grammar = r'''
!file: namelist*
!namelist: "&" NAME pair* "/"
!pair: NAME "=" /\S+/

NAME: /[a-zA-Z_0-9!]+/

%import common.WS
%ignore WS
'''

    def __init__(self):
        self.parser: Lark = Lark(
            grammar=NamelistGrammar.grammar,
            start='file',
            parser='lalr',
        )
        self.tree: Tree = None 
        self.options: dict = {}
        self.ran_block: bool = False

    def parse(self, text: str):
        self.tree = self.parser.parse(text=text)

    def value_filter(self, obj):
        # Looks convoluted. But not sure of a better solution.
        try:
            obj = int(obj)
        except ValueError:
            try:
                obj = float(obj)
            except ValueError:  # Must be string.
                obj = str(obj)
        
        return obj

    def get_namelists(self) -> dict:
        namelists = None

        for namelist in self.tree.find_data(data='namelist'):
            if namelists is None: namelists = {}

            dict_name = namelist.children[1].value.lower()
            dict_value = {}
            
            for pair in namelist.find_data(data='pair'):
                key = pair.children[0].value
                value = pair.children[2].value
                dict_value[key] = self.value_filter(value)

            namelists[dict_name] = dict_value

        return namelists

class QeFileGrammar:
    grammar = r'''
!file: (namelist | block)*
!namelist: "&" NAME pair* "/"
!pair: NAME "=" /\S+/
!block: block_title block_content
!block_title: "ATOMIC_SPECIES"
    | "CELL_PARAMETERS"
    | "ATOMIC_POSITIONS"
    | "K_POINTS"
    | "OCCUPATIONS"
!block_option: "angstrom"
    | "bohr"
    | "alat"
    | "crystal"
    | "crystal_b"
    | "automatic"
!block_content: /(.|\n)+?\n *\n/

NAME: /[a-zA-Z_0-9!]+/

%import common.WS
%ignore WS
'''

    def __init__(self):
        self.parser: Lark = Lark(
            grammar=QeFileGrammar.grammar,
            start='file',
            parser='lalr',
        )
        self.tree: Tree = None 
        self.options: dict = {}
        self.ran_block: bool = False

    def parse(self, text: str):
        self.tree = self.parser.parse(text=text)

    def value_filter(self, obj):
        # Looks convoluted. But not sure of a better solution.
        try:
            obj = int(obj)
        except ValueError:
            try:
                obj = float(obj)
            except ValueError:
                if obj=='.true.':
                    obj = True
                elif obj=='.false.':
                    obj = False
                else: # Must be a string.
                    obj = obj.removeprefix("'").removesuffix("'")
        
        return obj

    def get_namelists(self) -> dict:
        namelists = None

        for namelist in self.tree.find_data(data='namelist'):
            if namelists is None: namelists = {}

            dict_name = namelist.children[1].value.lower()
            dict_value = {}
            
            for pair in namelist.find_data(data='pair'):
                key = pair.children[0].value
                value = pair.children[2].value
                dict_value[key] = self.value_filter(value)

            namelists[dict_name] = dict_value

        return namelists

    def gen_list_of_lists(self, input_list: List[str]):
        output_list: list = []
        
        for list_item in input_list:
            filtered_list = [self.value_filter(item.strip()) for item in list_item.split(sep=' ') if item.strip() != '']
            if len(filtered_list)>0:
                if len(filtered_list)==6:   # case of K_POINTS automatic \n <> <> <> 0 0 0
                    output_list.append(filtered_list[:3])
                elif len(filtered_list)==0: # case of <num_kpoints>
                    continue
                else:
                    output_list.append(filtered_list)

        return output_list

    def get_blocks(self) -> dict:
        blocks = None
        for item in self.tree.find_data(data='block'):
            if blocks is None: blocks = {}

            block_key = item.children[0].children[0].value.lower()
            block_lines = item.children[1].children[0].value.splitlines()
            block_value = None

            # Add options.
            if block_key=='atomic_species':
                block_value = self.gen_list_of_lists(block_lines)
            elif block_key=='cell_parameters':
                self.options['cell_units'] = block_lines[0].strip()
                block_value = self.gen_list_of_lists(block_lines[1:])
            elif block_key=='atomic_positions':
                self.options['position_units'] = block_lines[0].strip()
                block_value = self.gen_list_of_lists(block_lines[1:])
            elif block_key=='k_points':
                self.options['kpoints_type'] = block_lines[0].strip()
                block_key = 'kpoints'
                block_value = self.gen_list_of_lists(block_lines[1:])
                # Accounting for cases like [[2, 2, 2]]
                if len(block_value)==1 and isinstance(block_value[0], list):
                    block_value = block_value[0]
                if isinstance(block_value[0], list) and len(block_value[0])==1:
                    del block_value[0]
                for row_idx, row in enumerate(block_value):
                    if isinstance(row, list) and isinstance(row[-1], str):
                        block_value[row_idx][-1] = row[-1].lstrip('!')

            # Add block entry.
            blocks[block_key] = block_value

        self.ran_block = True

        return blocks 

    def get_options(self) -> dict:
        if not self.ran_block: self.get_blocks()
        
        return self.options

class IbravType:
    SC = 1
    FCC = 2
    BCC = 3
    TETRA = 6
    ORTHO = 8
    HEX = 4
    FREE = 0

    def __init__(self, input_dict: dict):
        self.input_dict: dict = input_dict

        self.ibrav_str: str = 'free'
        if self.input_dict['atoms']['write_cell_units']=='alat':        
            self.ibrav_str: str = self.input_dict.get('atoms', {}).get('read', {}).get('cell', {}).get('alat_info', {}).get('ibrav')
            if self.ibrav_str is None: self.ibrav_str = 'free'

    def get_idx(self):
        if self.ibrav_str=='sc':
            return IbravType.SC
        elif self.ibrav_str=='fcc':
            return IbravType.FCC
        elif self.ibrav_str=='bcc':
            return IbravType.BCC
        elif self.ibrav_str=='tetra':
            return IbravType.TETRA
        elif self.ibrav_str=='ortho':
            return IbravType.ORTHO
        elif self.ibrav_str=='hex':
            return IbravType.HEX
        elif self.ibrav_str=='free':
            return IbravType.FREE
        else:
            return IbravType.FREE

class QePwInputFile:
    def __init__(self, qeinpdict: dict, input_dict: dict = None):
        self.qeinpdict: dict = qeinpdict
        self.input_dict: dict = input_dict

    def get_type_str(self, item):
        if isinstance(item, str):
            return f"'{item}'"
        elif isinstance(item, bool):
            return '.true.' if True else '.false'
        else:
            return item

    @staticmethod
    def read_pw(string_content: str) -> dict:
        pw_dict = {}
        parser = QeFileGrammar()
        parser.parse(text=string_content)

        if (blocks:= parser.get_blocks()) is not None:
            pw_dict['blocks'] = blocks

        if (blocks:= parser.get_namelists()) is not None:
            pw_dict['namelists'] = blocks

        pw_dict.update(parser.get_options())

        return pw_dict

    @staticmethod
    def read_general(string_content: str) -> dict:
        '''
        file consists of only namelists in this case.
        Applicable to ph.x, dos.x, epw.x and such.
        '''
        parser = NamelistGrammar()
        parser.parse(text=string_content)
        output_dict = {
            'namelists': parser.get_namelists()
        }

        return output_dict

    @staticmethod
    def write_general(general_dict: dict) -> str:
        output = ''
        
        # Write namelists.
        if general_dict.get('namelists') is not None:
            for key, value in general_dict['namelists'].items():
                output += f'&{key.upper()}\n'
                for value_key, value_value in value.items():
                    output += f'{value_key}={value_value}\n'
                output += '/\n'

        # Write blocks.
        if general_dict.get('blocks') is not None:
            for key, value in general_dict['blocks'].items():
                output += f'{key.upper()} '

                if key.lower() is 'cell_parameters':
                    output += f'{general_dict["cell_units"]}\n'
                
                if key.lower() is 'atomic_positions':
                    output += f'{general_dict["position_units"]}\n'

                if key.lower() is 'k_points':
                    output += f'{general_dict["kpoints_type"]}\n'

                for row in value:
                    for col in row:
                        output += f'{col} '
                    output += '\n'
                output += '\n'

        return output

    def get_control(self):
        output = ''
        output += '&CONTROL\n'

        if self.qeinpdict['namelists']['control'] is not None:
            for key, value in self.qeinpdict['namelists']['control'].items():
                output += f'{key}={self.get_type_str(value)}\n'

        output += '/\n\n'

        return output

    def get_system(self):
        # preprocess.
        if self.qeinpdict['cell_units']=='alat':
            ibrav_str = self.qeinpdict['namelists']['system']['ibrav']

            if ibrav_str=='sc':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
            elif ibrav_str=='fcc':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
            elif ibrav_str=='bcc':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
            elif ibrav_str=='tetra':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
                self.qeinpdict['namelists']['system']['C'] = self.input_dict['atoms']['read']['cell']['alat_info']['C']
            elif ibrav_str=='ortho':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
                self.qeinpdict['namelists']['system']['B'] = self.input_dict['atoms']['read']['cell']['alat_info']['B']
                self.qeinpdict['namelists']['system']['C'] = self.input_dict['atoms']['read']['cell']['alat_info']['C']
            elif ibrav_str=='hex':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
                self.qeinpdict['namelists']['system']['C'] = self.input_dict['atoms']['read']['cell']['alat_info']['C']
            elif ibrav_str=='free':
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']
            else:
                self.qeinpdict['namelists']['system']['A'] = self.input_dict['atoms']['read']['cell']['alat_info']['A']

        output = ''
        output += '&SYSTEM\n'

        if self.qeinpdict['namelists']['system'] is not None:
            for key, value in self.qeinpdict['namelists']['system'].items():
                output += f'{key}={self.get_type_str(value)}\n'

        output += '/\n\n'

        return output

    def get_electrons(self):
        output = ''
        output += '&ELECTRONS\n'

        if self.qeinpdict['namelists']['electrons'] is not None:
            for key, value in self.qeinpdict['namelists']['electrons'].items():
                output += f'{key}={self.get_type_str(value)}\n'

        output += '/\n\n'

        return output

    def get_ions(self):
        output = ''
        output += '&IONS\n'

        if self.qeinpdict['namelists']['ions'] is not None:
            for key, value in self.qeinpdict['namelists']['ions'].items():
                output += f'{key}={self.get_type_str(value)}\n'

        output += '/\n\n'

        return output

    def get_cell(self):
        output = ''

        output += '&CELL\n'

        if self.qeinpdict['namelists']['cell'] is not None:
            for key, value in self.qeinpdict['namelists']['cell'].items():
                output += f'{key}={self.get_type_str(value)}\n'

        output += '/\n\n'

        return output

    def get_occupations(self):
        output = 'OCCUPATIONS\n'

        if self.qeinpdict['blocks']['occupations'] is not None:
            for row in self.qeinpdict['blocks']['occupations']:
                output += f'{row:15.10f}\n'


        output += '\n'

        return output

    def get_atomic_species(self):
        output = ''
        output += 'ATOMIC_SPECIES\n'

        if self.qeinpdict['blocks']['atomic_species'] is not None:
            for row in self.qeinpdict['blocks']['atomic_species']:
                output += f'{row[0]} {row[1]} {row[2]}\n'
        output += '\n'

        return output

    def get_cell_parameters(self):
        output = ''

        # If units is alat, need not print anything. 
        if self.qeinpdict['cell_units'] == 'alat':
            return output

        cell_units = self.qeinpdict["cell_units"]
        output += f'CELL_PARAMETERS {cell_units}\n'
        if self.qeinpdict['blocks']['cell_parameters'] is not None:
            for row in self.qeinpdict['blocks']['cell_parameters']:
                output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n'
        output += '\n'

        return output

    def get_atomic_positions(self):
        output = ''
        pos_units = self.qeinpdict['position_units']
        output += f'ATOMIC_POSITIONS {pos_units}\n'

        if self.qeinpdict['blocks']['atomic_positions'] is not None:
            pos_list = self.qeinpdict['blocks']['atomic_positions']
            for row in pos_list:
                for col in row:
                    if isinstance(col, float):
                        output += f'{col:15.10f} '
                    else:
                        output += f'{col} '
                output += '\n'
        output += '\n'

        return output

    def get_kpoints(self):
        kpoints_type = self.qeinpdict['kpoints_type']
        kpoints = self.qeinpdict['blocks']['kpoints']
        output = f'K_POINTS {kpoints_type}\n'

        if kpoints_type=='automatic':
            output += f'{int(kpoints[0])} {int(kpoints[1])} {int(kpoints[2])} 0 0 0\n'
        elif kpoints_type=='crystal':
            output += f'{len(kpoints)}\n'
            for row in kpoints:
                factor = '1.0' if len(row)==3 else f'{row[3]:15.10f}'
                output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} {factor}\n'
        elif kpoints_type=='crystal_b':
            output += f'{len(kpoints)}\n'
            for row in kpoints:
                output += f'{row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f} {row[3]} !{row[4]}\n'

        output += '\n'

        return output 

    def get_input_str(self):
        '''
        TODO: Refactor this to not include many class method calls. 
        '''
        output = ''

        if 'control' in self.qeinpdict['namelists']: output += self.get_control()
        if 'system' in self.qeinpdict['namelists']: output += self.get_system()
        if 'electrons' in self.qeinpdict['namelists']: output += self.get_electrons()
        if 'ions' in self.qeinpdict['namelists']: output += self.get_ions()
        if 'cell' in self.qeinpdict['namelists']: output += self.get_cell()
        if 'occupations' in self.qeinpdict['blocks']: output += self.get_occupations()
        if 'atomic_species' in self.qeinpdict['blocks']: output += self.get_atomic_species()
        if 'atomic_positions' in self.qeinpdict['blocks']: output += self.get_atomic_positions()
        if 'cell_parameters' in self.qeinpdict['blocks']: output += self.get_cell_parameters()
        if 'kpoints' in self.qeinpdict['blocks']: output += self.get_kpoints()

        return output
#endregions