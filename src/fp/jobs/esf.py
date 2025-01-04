#region: Modules.
from fp.inputs.input_main import Input
from fp.io.strings import write_str_2_f
from fp.flows.run import run_and_wait_command
import os 
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class EsfJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        
        self.job_esfxctph = \
f'''#!/bin/bash
{self.scheduler.get_sched_header(self.job_info)}

echo "\nStarting eph switched calculation"
write_eph_h5.py ./tmp struct {self.input.xctph.num_epw_qpts} {self.input.xctph.num_epw_cond_bands} {self.input.xctph.num_epw_val_bands} --switch_nu_and_cart
cp eph.h5 eph_esf.h5
echo "Done eph switched calculation\n"

echo "\nStarting xctph esf calculation"
compute_xctph.py ./eph_esf.h5 ./xct.h5 {self.input.xctph.num_exciton_states}  --add_electron_part --add_hole_part 
mv xctph.h5 xctph_esf.h5
echo "Done xctph esf calculation\n"

# Print stuff if needed. 
echo "\nStaring printing"
print_xctph.py ./xctph_esf.h5 --switch_nu_and_cart_eVA
mv xctph.dat xctph_esf.dat
echo "Done printing\n"
'''

        self.jobs = [
            './job_esfxctph.sh',
        ]

    def create(self):
        write_str_2_f('job_esfxctph.sh', self.job_esfxctph)

    def run(self, total_time):
        for job in self.jobs:
            total_time = run_and_wait_command(job, self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'job_esfxctph.sh',

            'xct.h5',
            'eph*.h5',
            'eph*.dat',
            'xctph*.h5',
            'xctph*.dat',
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'job_esfxctph.sh',

            'xct.h5',
            'eph*.h5',
            'eph*.dat',
            'xctph*.h5',
            'xctph*.dat',
        ] 


        for inode in inodes:
            os.system(f'rm -rf ./{inode}')
#endregion
