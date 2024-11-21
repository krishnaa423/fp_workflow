#region: Modules.
from fp.inputs.input_main import *
from fp.io.strings import *
from fp.flows.run import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class XctPhJob:
    def __init__(
        self,
        input: Input,
    ):
        self.input: Input = input
        
        self.job_xctph = \
f'''#!/bin/bash
{self.input.scheduler.get_sched_header(self.input.xctph.job_desc)}

rm -rf xctph.out
touch xctph.out

exec &> xctph.out

echo "\\nStarting xct calculation"
write_xct_h5.py ./bseq_for_xctph/Q_\*/eigenvectors.h5
echo "Done xct calculation\\n"

echo "\\nStarting eph calculation"
write_eph_h5.py ./save struct {self.input.xctph.num_epw_qpts} {self.input.xctph.num_epw_cond_bands} {self.input.xctph.num_epw_val_bands} {self.input.scf.num_val_bands} {self.input.epw.job_desc.nk}
mv eph.h5 eph_xctph.h5
echo "Done eph calculation\\n"

echo "\\nStarting xctph elec-only calculation"
compute_xctph.py ./eph_xctph.h5 ./xct.h5 {self.input.xctph.num_exciton_states} --add_electron_part
mv xctph.h5 xctph_elec.h5 
echo "Done xctph elec-only calculation\\n"

echo "\\nStarting xctph hole-only calculation"
compute_xctph.py ./eph_xctph.h5 ./xct.h5 {self.input.xctph.num_exciton_states}  --add_hole_part 
mv xctph.h5 xctph_hole.h5
echo "Done xctph hole-only calculation\\n"

echo "\\nStarting xctph elec+hole calculation"
compute_xctph.py ./eph_xctph.h5 ./xct.h5 {self.input.xctph.num_exciton_states}  --add_electron_part --add_hole_part 
mv xctph.h5 xctph_elhole.h5
echo "Done xctph elec+hole calculation\\n"

# Print stuff if needed. 
echo "\\nStaring printing"
print_eph.py ./eph_xctph.h5
mv eph.dat eph_xctph.dat
print_xctph.py ./xctph_elhole.h5
mv xctph.dat xctph_elhole.dat
echo "Done printing\\n"
'''

        self.jobs = [
            'job_xctph.sh',
        ]

    def create(self):
        write_str_2_f('job_xctph.sh', self.job_xctph)

    def run(self, total_time):
        total_time = run_and_wait_command('./job_xctph.sh', self.input, total_time)

        return total_time

    def save(self, folder):
        inodes = [
            'job_xctph.sh',

            'xct.h5',
            'eph*.h5',
            'eph*.dat'
            'xctph*.h5',
            'xctph*.dat'
        ] 

        for inode in inodes:
            os.system(f'cp -r ./{inode} {folder}')

    def remove(self):
        inodes = [
            'job_xctph.sh',

            'xct.h5',
            'eph*.h5',
            'eph*.dat',
            'xctph*.h5',
            'xctph*.dat',
        ] 

        for inode in inodes:
            os.system(f'rm -rf ./{inode}')

#endregion
