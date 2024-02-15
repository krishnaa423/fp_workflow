from ase import Atoms
from ase.build import make_supercell
from ase.io import read, write
from ase.data import atomic_numbers
import numpy as np
import json 


symbols = [
    'Si',
    'Si'
]
alat = 5.43
cell = alat*np.array([
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
])
numbers = [atomic_numbers[symbol] for symbol in symbols]
scaled_positions = np.array([
    (0.00, 0.00, 0.00),
    (0.25, 0.25, 0.25),
])

struct = Atoms(
    numbers=numbers,
    cell=cell,
    scaled_positions=scaled_positions,
    pbc=[True, True, True],
)
sc = np.array([2, 2, 2])

# If loading CIF file. 
#struct = read('<name>.cif')

# Make supecell and rattle if needed. 
# struct = make_supercell(struct, np.diag(sc))
struct.rattle()

# Save. 
# write('struct.cif', struct)
# write('struct.xsf', struct)

# Create dict and write it. 
gen_input = {
    "struct_alat": alat,
    "struct_numbers": struct.get_atomic_numbers().tolist(),
    "struct_symbols": struct.get_chemical_symbols(),
    "struct_cell": struct.get_cell().tolist(),
    "struct_positions": struct.get_positions().tolist(),
}

with open('gen_input.json', 'w') as f: json.dump(gen_input, f)

# Optionally add to input.json too. 
with open('input.json', 'r') as f: input = json.load(f)

input['struct_alat'] = gen_input['struct_alat']
input['struct_numbers'] = gen_input['struct_numbers']
input['struct_symbols'] = gen_input['struct_symbols']
input['struct_cell'] = gen_input['struct_cell']
input['struct_positions'] = gen_input['struct_positions']
input['fold2bloch_sc'] = sc.tolist()


with open('input.json', 'w') as f: json.dump(input, f)
