#region modules
from fp.inputs.input_main import Input
from lark import Lark, Tree
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class BgwFileGrammar:
    grammar = r'''
!file: (flag | block)*
!block.2: "begin" /(.|\n)+?end/ 
!flag: /(?!begin).+/

%import common.WS
%ignore WS
'''

    def __init__(self):
        self.parser: Lark = Lark(
            grammar=BgwFileGrammar.grammar,
            start='file',
            parser='lalr',
        )
        self.tree: Tree = None 

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

    def get_list_of_lists(self, input_list):
        output_list = []
        
        for row in input_list:
            cols = [ self.value_filter(col) for col in row.split(sep=' ') if col!='']
            output_list.append(cols)

        return output_list

    def get_maps(self) -> dict:
        maps = {}

        for item in self.tree.find_data(data='flag'):
            filtered_line = [self.value_filter(col) for col in item.children[0].value.strip().split(sep=' ') if col != '']

            key = filtered_line[0]
            value = None
            if len(filtered_line)==1:
                value = ''
            else:
                value = filtered_line[1:]
                if len(value)==1: value = value[0]

            maps[key] = value

        return maps

    def get_blocks(self) -> dict:
        blocks = {}

        for item in self.tree.find_data(data='block'):
            lines = item.children[1].value.splitlines()
            key = lines[0].strip()
            value = self.get_list_of_lists(lines[1:-1])
            blocks[key] = value

        return blocks

class BgwInputFile:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

    @staticmethod
    def read_general(general_str: str) -> dict:
        parser = BgwFileGrammar()
        parser.parse(text=general_str)

        output_dict = {
            'maps': parser.get_maps(),
            'blocks': parser.get_blocks()
        }

        return output_dict


    @staticmethod
    def write_general(general_dict: dict) -> str:
        output = ''
        
        # Write maps.
        if general_dict.get('maps') is not None:
            for key, value in general_dict['maps'].items():
                # for now values are only lists or single data entries. (Cannot be dictionaries).
                if isinstance(value, list):
                    output += f'{key} '
                    for list_item in value:
                        output += f'{list_item} '
                else:
                    output += f'{key} {value}'
                
                output += '\n'

        # Write blocks.
        if general_dict.get('blocks') is not None:
            for key, value in general_dict['blocks'].items():
                if value is not None and len(value)!= 0:
                    output += f'begin {key} \n'
                    for row in value:
                        for col in row:
                            output += f'{col} '
                        output += '\n'
                    output += 'end\n'

        return output

    @staticmethod
    def update_dict(inputs_dict: dict, args_type: str, dest_dict: dict) -> dict:
        if inputs_dict is not None:
            if args_type=='override':
                dest_dict = inputs_dict.copy()
            if args_type=='extra':
                if 'maps' in inputs_dict:
                    if dest_dict.get('maps') is None: dest_dict['maps'] = {}
                    dest_dict['maps'].update(inputs_dict['maps'])
                if 'blocks' in inputs_dict:
                    if dest_dict.get('blocks') is None: dest_dict['blocks'] = {}
                    dest_dict['blocks'].update(inputs_dict['blocks'])

        return dest_dict
#endregions
