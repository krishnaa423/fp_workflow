from ase import Atoms
from ase.io import read, write
from ase.build import make_supercell
from ase.data import atomic_masses, atomic_numbers
import numpy as np
import json 
import os 
from io import StringIO

# Debugging. 
#os.chdir('/mnt/c/Users/User/Documents/Application_programming/Python/Samples/FirstScribbles')

files = {
    "scf": [
        "scf.in",
        "job_scf.sh",
    ],
    "wfn": [
        "wfn.in",
        "wfn_pw2bgw.in",
        "parabands.inp",
        "job_wfn.sh",
        "job_wfn_pw2bgw.sh",
        "job_parabands.sh",
    ],
    "wfnq": [
        "wfnq.in",
        "wfnq_pw2bgw.in",
        "job_wfnq.sh",
        "job_wfnq_pw2bgw.sh",
    ],
    "wfn_fi": [
        "wfn_fi.abi",
        "abi2bgw.in"
        "job_wfn_fi.sh",
        "job_abi2bgw.sh",
    ],
    "dfpt": [],     # Included in the above for now. 
    "elph": [],     # Included in the above for now. 
    "elband": [
        "elband.abi",
        "job_elband.sh",
    ],   
    "fold2bloch": [
        "job_fold2bloch.sh",
    ],
    "wannier": [
        "wannier.abi",
        "job_wannier.sh",
    ],
    "phband": [
        "phband.abi",
        "mrgddb.abi",
        "mrgddb.files",
        "anaddb.abi",
        "anaddb.files",
        "job_phband.sh",
        "job_mrgddb.sh",
        "job_anaddb.sh",
    ],
    "gwbse": [
        "gwbse.abi",
        "job_gwbse.sh",
    ],
    "eps": [
        "epsilon.inp",
        "job_epsilon.sh",
    ],
    "sig": [
        "bands.in",
        "sigma.inp",
        "inteqp.inp",
        "job_bands.in",
        "job_sigma.sh",
        "job_inteqp.sh",
    ],
    "abs": [
        "kernel.inp",
        "absorption.inp",
        "plotxct.inp",
        "job_kernel.sh",
        "job_absorption.sh",
        "job_plotxct.sh",
    ],
    "esf": [
        "esf.py",
        "job_esf.sh"
    ],
    "step": [
        "step.py",
        "job_step.sh"
    ],
    "update": [
        "update.py",
        "job_update.sh"
    ],
}

# region: Variables
input: dict = None 
# endregion

# region: Functions
def write_struct():
    struct = Atoms(
        numbers=input["struct_numbers"],
        cell=input["struct_cell"],
        positions=input["struct_positions"],
        pbc=[True, True, True],
    )
    write('struct.xsf', struct)
    write('struct.cif', struct)

def write_string(filename, string):
    with open(filename, 'w') as f: f.write(string)
    
def write_string_exec(filename, string):
    with open(filename, 'w') as f: f.write(string)
    os.system(f'chmod u+x {filename}')
    
def get_input():
    global input
    
    with open('input.json', 'r') as f: input = json.load(f)
    
def get_scf_cell():
    # In alat units.
    
    cell: np.ndarray = np.array(input['struct_cell'])/input['struct_alat']
    
    cell_string = StringIO()
    
    for row in cell:
        for item in row:
            cell_string.write(f'{item:10.15f} ')
        cell_string.write('\n')
    
    return cell_string.getvalue()

def get_scf_atomic_species():
    
    symbols = np.unique(np.array(input["struct_symbols"])).tolist()
    
    atomic_species_string = StringIO()
    
    for species in symbols:
        atomic_species_string.write(f'{species} {atomic_masses[atomic_numbers[species]]} {species}_ONCV_PBE_sr.upf\n')        # Or fr if using spin-orbit version.
    
    return atomic_species_string.getvalue()

def get_scf_atomic_positions():
    '''
    Get the atomic position string for SCF file.
    Units are in Angstroms. 
    '''
    
    atomic_symbols = input['struct_symbols']
    atomic_positions = input['struct_positions']
    
    atomic_positions_string = StringIO()
    
    for atm_idx, atomic_symbol in enumerate(atomic_symbols):
        atomic_positions_string.write(f'{atomic_symbol} {atomic_positions[atm_idx][0]:15.10f} {atomic_positions[atm_idx][1]:15.10f} {atomic_positions[atm_idx][2]:15.10f}\n')
    
    return atomic_positions_string.getvalue() 
    
def get_scf_natoms():
    return str(len(input['struct_symbols']))

def get_scf_ntypes():
    return str(np.unique(np.array(input['struct_symbols'])).size)

def get_abinit_znucl():
    symbols = np.unique(np.array(input["struct_symbols"])).tolist()
    
    return " ".join(list(map(lambda x: str(atomic_numbers[x]), symbols)))

def get_abinit_pseudos():
    symbols = np.unique(np.array(input["struct_symbols"])).tolist()
    
    return " ".join(list(map(lambda x: f"{x}.psp8, ", symbols)))[:-2]       # Removes the last comma. 

def get_abinit_typat():
    symbols = np.array(input["struct_symbols"]).tolist()
    numbers = list(map(lambda x: atomic_numbers[x], symbols))
    
    unique_symbols = np.unique(np.array(input["struct_symbols"])).tolist()
    unique_numbers = list(map(lambda x: atomic_numbers[x], unique_symbols))
    
    typat = [unique_numbers.index(number)+1 for number in numbers]
    typat = list(map(str, typat))
    
    
    return " ".join(typat)

def get_abinit_xred():
    '''
    In reduced units. 
    '''
    
    cell = np.array(input['struct_cell'])
    positions = np.array(input['struct_positions'])
    
    xred_string = StringIO()
    
    for index, position in enumerate(positions):
        xred = np.matmul(np.linalg.inv(cell.transpose()),position)
        xred_string.write(f'{xred[0]:15.10f} {xred[1]:15.10f} {xred[2]:15.10f}\n')
    
    return xred_string.getvalue()

def get_sched_exec(prefix):
    if len(input['sched_exec']) != 0:
        if len(input["sched_launch"]) == 0: 
            tasks = input[f'{prefix}_tasks']
            return f"{input['sched_exec']} -n {tasks} "
        else:
            return input['sched_exec']
    else:
        return ""

# endregion

# region: Create variables
get_input()
write_struct()
# endregion

# region: scf
scf_string = \
f'''&CONTROL
calculation='scf'
outdir='./tmp'
prefix='struct'
pseudo_dir='./qe_pseudos'
tprnfor=.true.
/

&SYSTEM
ibrav=0
ntyp={get_scf_ntypes()}
nat={get_scf_natoms()}
!noncolin=.true.
!lspinorb=.true.
ecutwfc={input["scf_ecutwfc"]}
A={input["struct_alat"]}
/

&ELECTRONS
/

CELL_PARAMETERS alat
{get_scf_cell()}

ATOMIC_SPECIES
{get_scf_atomic_species()}

ATOMIC_POSITIONS angstrom
{get_scf_atomic_positions()}

K_POINTS automatic
1 1 1 0 0 0
'''
job_scf_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account=m3571
#{input["sched_bash_prefix"]} --constraint=cpu
#{input["sched_bash_prefix"]} --qos=debug
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["scf_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["scf_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["scf_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('scf_job')}pw.x < scf.in &> scf.out 
'''

write_string('scf.in', scf_string)
write_string_exec('job_scf.sh', job_scf_string)
# endregion

# region: wfn
wfn_string = \
f'''&CONTROL
calculation='bands'
outdir='./tmp'
prefix='struct'
pseudo_dir='./qe_pseudos'
tprnfor=.true.
/

&SYSTEM
ibrav=0
ntyp={get_scf_ntypes()}
nat={get_scf_natoms()}
nbnd={input["wfn_bands"]}
!noncolin=.true.
!lspinorb=.true.
ecutwfc={input["scf_ecutwfc"]}
A={input["struct_alat"]}
/

&ELECTRONS
/

CELL_PARAMETERS alat
{get_scf_cell()}

ATOMIC_SPECIES
{get_scf_atomic_species()}

ATOMIC_POSITIONS angstrom
{get_scf_atomic_positions()}

K_POINTS crystal
1
0.0 0.0 0.0 1.0
'''
wfn_pw2bgw_string = \
f'''&input_pw2bgw
outdir = './tmp'
prefix = 'struct'
real_or_complex = 2
wfng_flag = .true.
wfng_file = 'WFN'
wfng_kgrid = .true.
wfng_nk1 = 1
wfng_nk2 = 1
wfng_nk3 = 1
wfng_dk1 = 0.0
wfng_dk2 = 0.0
wfng_dk3 = 0.0
rhog_flag = .true.
rhog_file = 'RHO'
vxcg_flag = .true.
vxcg_file = 'VXC'
vxc_flag = .true.
vxc_file = 'vxc.dat'
vxc_diag_nmin = 1
vxc_diag_nmax = {input["wfn_bands"]}
vxc_offdiag_nmin = 0
vxc_offdiag_nmax = 0
vscg_flag = .true.
vscg_file = 'VSC'
vkbg_flag = .true.
vkbg_file = 'VKB'
/
'''
wfn_parabands_string = \
f'''input_wfn_file WFN
output_wfn_file WFN_generated
vkb_file VKB
vsc_file VSC

number_bands {input["wfn_parabands"]}
'''
job_wfn_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["wfn_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["wfn_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["wfn_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('wfn_job')}pw.x < wfn.in &> wfn.out 
'''
job_wfn_pw2bgw_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes=1
#{input["sched_bash_prefix"]} --ntasks=1
#{input["sched_bash_prefix"]} --time=00:30:00
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{input["sched_exec"]}pw2bgw.x -pd .true. < wfn_pw2bgw.in &> wfn_pw2bgw.out 
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB
# wfn2hdf.x BIN WFN WFN.h5
'''
job_parabands_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["wfn_job_parabands_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["wfn_job_parabands_tasks"]}
#{input["sched_bash_prefix"]} --time={input["wfn_job_parabands_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('wfn_job_parabands')}parabands.cplx.x &> parabands.out
mv WFN_generated WFN.h5
hdf2wfn.x BIN WFN.h5 WFN
'''

write_string('wfn.in', wfn_string)
write_string('wfn_pw2bgw.in', wfn_pw2bgw_string)
write_string('parabands.inp', wfn_parabands_string)
write_string_exec('job_wfn.sh', job_wfn_string)
write_string_exec('job_wfn_pw2bgw.sh', job_wfn_pw2bgw_string)
write_string_exec('job_parabands.sh', job_parabands_string)
# endregion

# region: wfnq
wfnq_string = \
f'''&CONTROL
calculation='bands'
outdir='./tmp'
prefix='struct'
pseudo_dir='./qe_pseudos'
tprnfor=.true.
/

&SYSTEM
ibrav=0
ntyp={get_scf_ntypes()}
nat={get_scf_natoms()}
nbnd={input["wfnq_bands"]}
!noncolin=.true.
!lspinorb=.true.
ecutwfc={input["scf_ecutwfc"]}
A={input["struct_alat"]}
/

&ELECTRONS
/

CELL_PARAMETERS alat
{get_scf_cell()}

ATOMIC_SPECIES
{get_scf_atomic_species()}

ATOMIC_POSITIONS angstrom
{get_scf_atomic_positions()}

K_POINTS crystal
1
0.0 0.0 0.001 1.0
'''
wfnq_pw2bgw_string = \
f'''&input_pw2bgw
outdir = './tmp'
prefix = 'struct'
real_or_complex = 2
wfng_flag = .true.
wfng_file = 'WFNq'
wfng_kgrid = .true.
wfng_nk1 = 1
wfng_nk2 = 1
wfng_nk3 = 1
wfng_dk1 = 0.0
wfng_dk2 = 0.0
wfng_dk3 = 0.001
/
'''
job_wfnq_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["wfnq_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["wfnq_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["wfnq_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('wfnq_job')}pw.x < wfnq.in &> wfnq.out 
'''
job_wfnq_pw2bgw_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes=1
#{input["sched_bash_prefix"]} --ntasks=1
#{input["sched_bash_prefix"]} --time=00:30:00
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{input["sched_exec"]}pw2bgw.x -pd .true. < wfnq_pw2bgw.in &> wfnq_pw2bgw.out 
cp ./tmp/WFNq ./WFNq
wfn2hdf.x BIN WFNq WFNq.h5
'''

write_string('wfnq.in', wfnq_string)
write_string('wfnq_pw2bgw.in', wfnq_pw2bgw_string)
write_string_exec('job_wfnq.sh', job_wfnq_string)
write_string_exec('job_wfnq_pw2bgw.sh', job_wfnq_pw2bgw_string)
# endregion

# region: wfn_fi
wfn_fi_string = \
f'''########################################################
# Top-level global control. 
########################################################

# IO. 
iomode3 3
prtwf 1
prtwf_full 1

# Datasets. 
ndtset 3
#jdtset             # Control what dataset calcs to run. 

# Symmetry. 
nsym 1

########################################################
# Structure. 
########################################################

# Unit cell vectors. 
acell 3*{input["struct_alat"]} Angstrom 
rprim
{get_scf_cell()} 

# Atom types. 
ntypat {get_scf_ntypes()}
znucl {get_abinit_znucl()}
pp_dirpath "./abinit_pseudos"
pseudos "{get_abinit_pseudos()}"

# Atoms. 
natom {get_scf_natoms()}
typat {get_abinit_typat()}
xred
{get_abinit_xred()}

# k-points. 
kptopt 3    # Generate k-points using ngkpt and kptrlatt without taking symmetry into account. 
ngkpt 1 1 1
nshiftk 1
shiftk
0.0 0.0 0.0
prtkpt 1

########################################################
# SCF. Dataset: 1
########################################################

# Bands and cutoff. 
#nband1 4     # Only valence bands. Comment out to automatically choose for SCF.  
ecut1 {input["wfn_fi_ecut"]}   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria. 
toldfe1 1.0d-12   # Stops when difference in energy is less that this amount. 

# Printing options. 
prtden1 1
prtvxc1 1

########################################################
# NSCF. Dataset: 2
########################################################

# Calc type and previous inputs.
iscf2 -2 
getden2 -1

# Band and ecut. 
nband2 {input["wfn_fi_bands"]}
ecut2 {input["wfn_fi_ecut"]}   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria.
tolwfr2  1.0d-12

########################################################
# DFPT and ELPH. Dataset: 3
########################################################

# Phonon calculation parameters and previous inputs. 
rfphon3 1       # Activate response calc w.r.t atomic displacements. 
getwfk3 -1      # Read wavefunction from NSCF calculation.  

# q-points. 
nqpt3 1
qpt3 0.0 0.0 0.0

# Bands and cutoff. 
nband3 {input["wfn_fi_bands"]}    # Only valence bands. 
ecut3 {input["wfn_fi_ecut"]}   # 12.0 Ha, which is 24.0 Ry. 

# Stopping criteria.
tolvrs3 1.0d-8       # Stopping criteria for self-consistent DFPT calc.  

# Printing options. 
#prtfull1wf3 1
prtgkk3 1

# Parallelization control. 
#paral_rf3 1
#nppert 40          # If there are 240 perturbations, then This is 240/40=6 perts per processor. 
'''
abi2bgw_string = \
f'''wfng_file_abi 'wfn_fio_DS2_FULL_WFK'
wfng_flag .true.
wfng_file 'WFN_fi'
wfng_nk1 1
wfng_nk2 1
wfng_nk3 1
wfng_dk1 0
wfng_dk2 0
wfng_dk3 0
rhog_file_abi 'wfn_fio_DS1_DEN'
rhog_flag .false.
rhog_file 'RHO'
cell_symmetry 0
symrel_file_flag .false.
vxcg_file_abi 'wfn_fio_DS1_VXC'
vxcg_flag .false.
vxcg_file 'VXC'
'''
job_wfn_fi_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["wfn_fi_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["wfn_fi_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["wfn_fi_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

rm -rf wfn_fi.abo* wfn_fio* wfn_fi.out
{get_sched_exec('wfn_fi_job')}abinit wfn_fi.abi &> wfn_fi.out 
'''
job_abi2bgw_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account=m3571
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes=1
#{input["sched_bash_prefix"]} --ntasks=1
#{input["sched_bash_prefix"]} --time=00:30:00
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{input["sched_exec"]}abi2bgw.x &> abi2bgw.out
wfn2hdf.x BIN WFN_fi WFN_fi.h5
'''

write_string('wfn_fi.abi', wfn_fi_string)
write_string('abi2bgw.inp', abi2bgw_string)
write_string_exec('job_wfn_fi.sh', job_wfn_fi_string)
write_string_exec('job_abi2bgw.sh', job_abi2bgw_string)
# endregion

# region: elband
elband_string = \
f'''########################################################
# Top-level global control. 
########################################################

# IO. 
iomode 3

# Datasets. 
ndtset 2
#jdtset             # Control what dataset calcs to run. 

# Symmetry. 
#nsym 1

########################################################
# Structure. 
########################################################

# Unit cell vectors. 
acell 3*{input["struct_alat"]}
rprim
{get_scf_cell()}

# Atom types. 
ntypat {get_scf_ntypes()}
znucl {get_abinit_znucl()}
pp_dirpath "./abinit_pseudos"
pseudos "{get_abinit_pseudos()}"

# Atoms. 
natom {get_scf_natoms()}
typat {get_abinit_typat()}
xred
{get_abinit_xred()}



########################################################
# SCF. Dataset: 1
########################################################

# k-points. 
kptopt1 3    # Generate k-points using ngkpt and kptrlatt without taking symmetry into account. 
ngkpt1 {input["elband_kpts"]}
nshiftk1 1
shiftk1
0.0 0.0 0.0
prtkpt1 1

# Bands and cutoff. 
#nband1 {input["wfn_fi_bands"]}     # Only valence bands. 
ecut1 {input["wfn_fi_ecut"]}   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria. 
toldfe1 1.0d-6   # Stops when difference in energy is less that this amount. 

# Printing options. 
prtden1 1
prtvxc1 1

########################################################
# NSCF. Dataset: 2
########################################################

# Calc type and previous inputs.
iscf2 -2 
getden2 -1

# k-points. 
kptopt2 -3
ndivsm2 10
kptbounds2
0.500 0.500 0.500 #L
0.000 0.000 0.000 #G
0.500 0.000 0.500 #X
0.375 0.375 0.750 #K

# Band and ecut. 
nband2 {input["wfn_fi_bands"]}
ecut2 {input["wfn_fi_ecut"]}   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria.
tolwfr2  1.0d-12
'''
job_elband_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["elband_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["elband_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["elband_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

rm -rf elband.abo* elbando* elband.out
{get_sched_exec('elband_job')}abinit elband.abi &> elband.out 
'''

write_string('elband.abi', elband_string)
write_string_exec('job_elband.sh', job_elband_string)
# endregion

# region: fold2Bloch
job_fold2bloch_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["fold2bloch_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["fold2bloch_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["fold2bloch_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('fold2bloch_job')}fold2Bloch elbando_DS2_WFK.nc {input["struct_supercell"][0]}:{input["struct_supercell"][1]}:{input["struct_supercell"][2]}
'''

write_string_exec('job_fold2bloch.sh', job_fold2bloch_string)
# endregion

# region: wannier
wannier_string = \
f'''########################################################
# Top-level global control. 
########################################################

# IO. 
iomode 3

# Datasets. 
ndtset 2
#jdtset             # Control what dataset calcs to run. 

# Symmetry. 
#nsym 1

########################################################
# Structure. 
########################################################

# Unit cell vectors. 
acell 3*10.18
rprim
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0

# Atom types. 
ntypat 1
znucl 14
pp_dirpath "./abinit_pseudos"
pseudos "Si.psp8"

# Atoms. 
natom 2
typat 1 1
xred
0.00 0.00 0.00
0.25 0.25 0.25



########################################################
# SCF. Dataset: 1
########################################################

# k-points. 
kptopt1 3    # Generate k-points using ngkpt and kptrlatt without taking symmetry into account. 
ngkpt1 2 2 2
nshiftk1 1
shiftk1
0.0 0.0 0.0
prtkpt1 1

# Bands and cutoff. 
nband1 4     # Only valence bands. 
ecut1 12.0   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria. 
toldfe1 1.0d-6   # Stops when difference in energy is less that this amount. 

# Printing options. 
prtden1 1
prtvxc1 1

########################################################
# NSCF and Wannier. Dataset: 2
########################################################

# Calc type and previous inputs.
iscf2 -2 
getden2 -1

# k-points. 
kptopt2 3    # Generate k-points using ngkpt and kptrlatt without taking symmetry into account. 
ngkpt2 2 2 2
nshiftk2 1
shiftk2
0.0 0.0 0.0
prtkpt2 1

# Band and ecut. 
nband2 10
ecut2 12.0   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria.
tolwfr2  1.0d-12

# Wannier. 
prtwant2 1
w90prtunk2 1
'''
job_wannier_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["wannier_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["wannier_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["wannier_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

rm -rf wannier.abo* wanniero* wannier.out
{get_sched_exec('wannier_job')}abinit wannier.abi &> wannier.out
'''

write_string('wannier.abi', wannier_string)
write_string_exec('job_wannier.sh', job_wannier_string)
# endregion

# region: phband
phband_string = \
f'''########################################################
# Top-level global control. 
########################################################

# IO. 
iomode 3

# Datasets. 
ndtset 3
#jdtset             # Control what dataset calcs to run. 

# Symmetry. 
#nsym 1

########################################################
# Structure. 
########################################################

# Unit cell vectors. 
acell 3*10.18
rprim
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0

# Atom types. 
ntypat 1
znucl 14
pp_dirpath "./abinit_pseudos"
pseudos "Si.psp8"

# Atoms. 
natom 2
typat 1 1
xred
0.00 0.00 0.00
0.25 0.25 0.25

# k-points. 
kptopt 3    # Generate k-points using ngkpt and kptrlatt without taking symmetry into account. 
ngkpt 1 1 1
nshiftk 1
shiftk
0.0 0.0 0.0
prtkpt 1

########################################################
# SCF. Dataset: 1
########################################################

# Bands and cutoff. 
nband1 4     # Only valence bands. 
ecut1 12.0   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria. 
toldfe1 1.0d-6   # Stops when difference in energy is less that this amount. 

# Printing options. 
prtden1 1
prtvxc1 1

########################################################
# NSCF. Dataset: 2
########################################################

# Calc type and previous inputs.
iscf2 -2 
getden2 -1

# Band and ecut. 
nband2 10
ecut2 12.0   # 12.0 Ha, which is 24.0 Ry. 

# Stop criteria.
tolwfr2  1.0d-12

########################################################
# DFPT and ELPH. Dataset: 3
########################################################

# Phonon calculation parameters and previous inputs. 
rfphon3 1       # Activate response calc w.r.t atomic displacements. 
getwfk3 -1      # Read wavefunction from NSCF calculation.  

# q-points. 
nqpt3 1
qpt3 0.0 0.0 0.0

# Bands and cutoff. 
nband3 10     # Only valence bands. 
ecut3 12.0   # 12.0 Ha, which is 24.0 Ry. 

# Stopping criteria.
tolvrs3 1.0d-8       # Stopping criteria for self-consistent DFPT calc.  

# Printing options. 
#prtfull1wf 1
prtgkk3 1

# Parallelization control. 
#paral_rf3 1
#nppert 40          # If there are 240 perturbations, then This is 240/40=6 perts per processor. 
'''
mrgddb_input_string = \
f'''asd
'''
mrgddb_files_string = \
f'''asd
'''
anaddb_input_string = \
f'''asd
'''
anaddb_files_string = \
f'''asd
'''
job_phband_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["phband_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["phband_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["phband_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

rm -rf phband.abo* phbando* phband.out
{get_sched_exec('phband_job')}abinit phband.abi &> phband.out 
'''
job_mrgddb_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes=1
#{input["sched_bash_prefix"]} --ntasks=1
#{input["sched_bash_prefix"]} --time=00:30:00
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{input["sched_exec"]}mrgddb < mrgddb.files &> mrgddb.out  
'''
job_anaddb_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["phband_job_anaddb_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["phband_job_anaddb_tasks"]}
#{input["sched_bash_prefix"]} --time={input["phband_job_anaddb_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('phband_job_anaddb')}anaddb < anaddb.files &> anaddb.out  
'''

write_string('phband.abi', phband_string)
write_string('mrgddb.abi', mrgddb_input_string)
write_string('mrgddb.files', mrgddb_files_string)
write_string('anaddb.abi', anaddb_input_string)
write_string('anaddb.files', anaddb_files_string)
write_string_exec('job_phband.sh', job_phband_string)
write_string_exec('job_mrgddb.sh', job_mrgddb_string)
write_string_exec('job_anaddb.sh', job_anaddb_string)
# endregion

# region: gwbse
gwbse_string = \
f'''as
'''
job_gwbse_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes=1
#{input["sched_bash_prefix"]} --ntasks=1
#{input["sched_bash_prefix"]} --time=00:30:00
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{input["sched_exec"]}abinit gwbse.abi &> gwbse.out 
'''

write_string('gwbse.abi', gwbse_string)
write_string_exec('job_gwbse.sh', job_gwbse_string)
# endregion

# region: eps
eps_string = \
f'''epsilon_cutoff {input["eps_cutoff"]}
number_bands {input["eps_bands"]}
use_wfn_hdf5
degeneracy_check_override

begin qpoints
0.0 0.0 0.001 1.0 1
end
'''
job_eps_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["eps_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["eps_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["eps_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('eps_job')}epsilon.cplx.x &> epsilon.out
'''

write_string('epsilon.inp', eps_string)
write_string_exec('job_epsilon.sh', job_eps_string)
# endregion

# region: sig
sig_string = \
f'''band_index_min {input["sig_bands"][0]}
band_index_max {input["sig_bands"][1]}
number_bands {input["eps_bands"]}

use_wfn_hdf5
degeneracy_check_override
no_symmetries_q_grid

begin kpoints
0.0 0.0 0.0 1.0
end
'''
job_sig_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["sig_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["sig_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["sig_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

ln -sf WFN.h5 WFN_inner.h5
{get_sched_exec('sig_job')}sigma.cplx.x &> sigma.out
'''

write_string('sigma.inp', sig_string)
write_string_exec('job_sigma.sh', job_sig_string)
# endregion

# region: abs
kernel_string = \
f'''number_val_bands {input["abs_coarse_bands"][0]}
number_cond_bands {input["abs_coarse_bands"][1]}

use_wfn_hdf5

no_symmetries_coarse_grid

screening_semiconductor
'''
absorption_string = \
f'''diagonalization

use_wfn_hdf5
degeneracy_check_override

number_val_bands_coarse {input["abs_coarse_bands"][0]}
number_val_bands_fine {input["abs_fine_bands"][0]}
number_cond_bands_coarse {input["abs_coarse_bands"][1]}
number_cond_bands_fine {input["abs_fine_bands"][1]}

screening_semiconductor

# use_velocity

use_momentum
polarization 0.0 0.0 0.001

gaussian_broadening
energy_resolution 0.15

dump_bse_hamiltonian

no_symmetries_coarse_grid
no_symmetries_fine_grid
no_symmetries_shifted_grid

eqp_co_corrections
write_eigenvectors 10
'''
plotxct_string = \
f'''hole_position {input["plotxct_hole_position"][0]} {input["plotxct_hole_position"][1]} {input["plotxct_hole_position"][2]} 
#plot_spin 1
plot_state {input["plotxct_state"]}
q_shift 0.0 0.0 0.001
supercell_size 1 1 1 
verbosity 2
use_wfn_hdf5

#spinor
#hole_spin 1
#electron_spin 2
'''
job_kernel_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["abs_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["abs_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["abs_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

ln -sf WFN.h5 WFN_co.h5
{get_sched_exec('abs_job')}kernel.cplx.x &> kernel.out
'''
job_absorption_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["abs_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["abs_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["abs_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('abs_job')}absorption.cplx.x &> absorption.out
'''
job_plotxct_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["abs_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["abs_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["abs_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('abs_job')}plotxct.cplx.x &> plotxct.out
volume.py scf.in espresso *.a3Dr a3dr xct_1_0.0_0.0_0.0.xsf xsf false abs2 true
rm *a3Dr
'''

write_string('kernel.inp', kernel_string)
write_string('absorption.inp', absorption_string)
write_string('plotxct.inp', plotxct_string)
write_string_exec('job_kernel.sh', job_kernel_string)
write_string_exec('job_absorption.sh', job_absorption_string)
write_string_exec('job_plotxct.sh', job_plotxct_string)
# endregion

# region: esf
job_esf_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["esf_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["esf_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["esf_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('esf_job')}python3 esf.py &> esf.out 
'''

write_string_exec('job_esf.sh', job_esf_string)
# endregion

# region: step
job_step_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["step_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["step_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["step_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('step_job')}python3 step.py &> step.out 
'''

write_string_exec('job_step.sh', job_step_string)
# endregion

# region: update
job_update_string = \
f'''#!/bin/bash
#{input["sched_bash_prefix"]} --account={input["sched_account"]}
#{input["sched_bash_prefix"]} --constraint={input["sched_constraint"]}
#{input["sched_bash_prefix"]} --qos={input["sched_qos"]}
#{input["sched_bash_prefix"]} --job-name=struct
#{input["sched_bash_prefix"]} --nodes={input["update_job_nodes"]}
#{input["sched_bash_prefix"]} --ntasks={input["update_job_tasks"]}
#{input["sched_bash_prefix"]} --time={input["update_job_time"]}
#{input["sched_bash_prefix"]} --mail-user=krishnaa.vadivel@yale.edu
#{input["sched_bash_prefix"]} --mail-type=ALL

{get_sched_exec('update_job')}python3 update.py &> update.out 
'''

write_string_exec('job_update.sh', job_update_string)
# endregion