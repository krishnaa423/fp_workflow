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
            match key:
                case 'exclude_bands':
                    # Skip if the value is None.
                    if value is None:
                        continue

                    # Continue is the value is not None.
                    output += f'{key} = '
                    last_idx = len(value)-1
                    for list_idx, list_item in enumerate(value):
                        if isinstance(list_item, tuple):
                            output += f' {list_item[0]}-{list_item[1]} ' if list_idx==last_idx else f' {list_item[0]}-{list_item[1]}, '
                        else:
                            output += f' {list_item} ' if list_idx==last_idx else f' {value}, '
                    output += '\n'
                case _ :      
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
