#region modules
from fp.inputs.input_main import Input
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class BgwInputFile:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input

    @staticmethod
    def read_general(general_str: str) -> dict:
        pass

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
