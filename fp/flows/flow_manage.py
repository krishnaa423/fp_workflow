#region: Modules.
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class FlowManage:
    def __init__(
        self,
        list_of_steps,
    ):
        self.list_of_steps: list = list_of_steps

    def create(self):
        for step in self.list_of_steps:
            step.create()

    def run(self, total_time=0.0):
        total_time: float = total_time

        for step in self.list_of_steps:
            total_time = step.run(total_time)

        # Write the total workflow run time. 
        print(f'Done whole worflow in {total_time:15.10f} seconds.\n\n', flush=True)

    def save(self, folder):       
        for step in self.list_of_steps:
            step.save(folder)

    def remove(self):
        for step in self.list_of_steps:
            step.remove()
         
#endregion