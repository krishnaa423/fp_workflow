#region: Modules.
import fp.schedulers as schedulers
from fp.io.strings import write_str_2_f
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class JobProcDesc:
    def __init__(
        self,
        nodes: int = None,
        ntasks: int = None,
        time: int = None,
        nk: int = None,
        ni: int = None,
        **kwargs,
    ):
        self.nodes: int = nodes
        self.ntasks: int = ntasks
        self.time: str = time
        self.nk: int = nk
        self.ni: int = ni

        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def from_job_id(id: str, input_dict: dict, ):
        assert id in input_dict['job_types'], f'id: {id} is not in job_types'

        return JobProcDesc(**input_dict['job_types'][id])


class Scheduler:
    def __init__(
        self,
        sched_dict: dict,
        **kwargs,
    ):
        self.sched_dict: dict = sched_dict

        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def from_input_dict(input_dict: dict, runtype='mainrun'):
        # sched_name = list(input_dict['scheduler'][runtype].keys())[0]
        # sched_class = getattr(schedulers, sched_name)
        
        # return sched_class(**input_dict['scheduler'][runtype][sched_name])

        key = input_dict['scheduler'][runtype]

        return Scheduler(sched_dict=input_dict['scheduler'][key])

    def add_header_commands(self):
        output = ''
        if self.sched_dict.get('header_commands') is not None:
            output += '\n'
            output += f'{self.sched_dict['header_commands']}'
            output += '\n'

        return output 
    
    def get_sched_header(self, job_desc: JobProcDesc):
        header = ''

        # Return if it is in interactive mode.
        if self.sched_dict.get('options', {}).get('is_interactive') is not None:
            if self.sched_dict['options']['is_interactive']:
                header += self.add_header_commands()
                return header
            
        if self.sched_dict.get('header') is not None:
            for key, value in self.sched_dict['header'].items():
                header += f"#{self.sched_dict['launch'].upper()} --{key}={value}\n"

        header += f"#{self.sched_dict['launch'].upper()} --nodes={job_desc.nodes}\n"
        header += f"#{self.sched_dict['launch'].upper()} --time={job_desc.time}\n"

        header += self.add_header_commands()

        return header

    def get_sched_mpi_prefix(self, job_desc: JobProcDesc):
        prefix = ''

        # If there is no mpi, just return empty string.
        if self.sched_dict.get('mpi') is None:
            return prefix
        if self.sched_dict['mpi']=='':
            return prefix

        prefix += f"{self.sched_dict['mpi']} "
        prefix += f'-n {job_desc.ntasks} '
        
        if self.sched_dict.get('options', {}).get('is_gpu') is not None:
            if self.sched_dict['options']['is_gpu']:
                prefix += f"--gpus-per-task={self.sched_dict['node_info']['gpus']} "

        return prefix
    
    def get_sched_mpi_infix(self, job_desc: JobProcDesc):
        infix = ''

        # If there is no mpi, just return empty string.
        if self.sched_dict.get('mpi') is None:
            return infix
        if self.sched_dict['mpi']=='':
            return infix

        if job_desc.nk is not None:
            infix += f' -nk {job_desc.nk}'

        if job_desc.ni is not None:
            infix += f' -ni {job_desc.ni}'

        return infix

    def get_sched_submit(self):
        submit = ''

        if self.sched_dict.get('launch') is not None:
            submit += f"{self.sched_dict['launch']} "

        return submit
    
    def create_interactive(self):
        file_contents = '#!/bin/bash\n'

        file_contents += f"{self.sched_dict['launch']} "

        if self.sched_dict.get('interactive', {}).get('args') is not None:
            for key, value in self.sched_dict['interactive']['args'].items():
                if value is None:
                    file_contents += f' --{key} '
                else:
                    file_contents += f' --{key}=={value} '


        if self.sched_dict.get('interactive', {}).get('extra_string') is not None:
            file_contents += f" {self.sched_dict['interactive']['extra_string']}"

        file_contents += '\n'

        write_str_2_f('job_interactive.sh', file_contents)
#endregion
