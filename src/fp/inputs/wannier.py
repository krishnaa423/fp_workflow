#region: Modules.
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class WannierWinFile:
    def __init__(self):
        pass

    @staticmethod
    def write_general(wan_dict: dict) -> str:
        output = ''
        
        # maps.
        for key, value in wan_dict['maps'].items():
            if isinstance(value, list):
                output += f'{key} = '
                for item in value:
                    output += f'{item} '
                output += '\n'
            else:
                output += f'{key} = {value}\n'

        # blocks.
        for key, value in wan_dict['blocks'].items():
            if value is not None:
                output += f'begin {key}\n'
                for row in value:
                    for col in row:
                        output += f'{col} '
                    output += '\n'
                output += f'end {key}\n\n'

        return output 

    @staticmethod
    def read_general(wan_str: str) -> dict:
        # TODO: if we want to read it.
        pass

#endregion
