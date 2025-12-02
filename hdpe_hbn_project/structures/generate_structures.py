import os
import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import add_vacuum

num_configs = 20        
num_monomers = 6        
slab_a = 2.50           
slab_c = 6.66           
slab_size = 4           
vacuum = 15.0           
output_dir = 'hdpe_hbn_project/structures'
os.makedirs(output_dir, exist_ok=True)

a1 = np.array([slab_a, 0, 0])
a2 = np.array([slab_a/2, slab_a*np.sqrt(3)/2, 0])
positions = [np.array([0.0, 0.0, 0.0]), np.array([slab_a/3, slab_a/3, 0.0])]  
symbols = ['B', 'N']

cell = np.array([
    [a1[0]*slab_size, a1[1]*slab_size, 0.0],
    [a2[0]*slab_size, a2[1]*slab_size, 0.0],
    [0.0, 0.0, slab_c]
])

atoms = Atoms(symbols=symbols * slab_size**2,
              positions=[pos + i*a1 + j*a2
                         for i in range(slab_size)
                         for j in range(slab_size)
                         for pos in positions],
              cell=cell,
              pbc=[True, True, False])

add_vacuum(atoms, vacuum)

cc_bond = 1.54  
ch_bond = 1.09  #
angle = 109.5 * np.pi / 180.0

hdpe_chain = Atoms()
for i in range(num_monomers):
    cx = i * cc_bond
    cy = 0.0
    cz = 0.0
    # hydrogens- tetrahedral
    h1 = [0.0, ch_bond * np.cos(angle/2), ch_bond * np.sin(angle/2)]
    h2 = [0.0, -ch_bond * np.cos(angle/2), -ch_bond * np.sin(angle/2)]
    hdpe_chain += Atoms(['C','H','H'],
                        positions=[[cx, cy, cz],
                                   [cx+h1[0], cy+h1[1], cz+h1[2]],
                                   [cx+h2[0], cy+h2[1], cz+h2[2]]])

com = hdpe_chain.get_center_of_mass()
hdpe_chain.translate(-np.array([com[0], com[1], 0.0]))
chain_length = (num_monomers - 1) * cc_bond

cell_x = np.linalg.norm(atoms.get_cell()[0])
cell_y = np.linalg.norm(atoms.get_cell()[1])
slab_z_max = atoms.positions[:,2].max()

for i in range(num_configs):
    slab_copy = atoms.copy()
    hdpe_copy = hdpe_chain.copy()

    # random small rotations 
    hdpe_copy.translate(-hdpe_copy.get_center_of_mass())
    hdpe_copy.rotate(np.random.rand()*360, 'x', rotate_cell=False)
    hdpe_copy.rotate(np.random.rand()*360, 'y', rotate_cell=False)
    hdpe_copy.rotate(np.random.rand()*360, 'z', rotate_cell=False)
    hdpe_copy.translate(hdpe_copy.get_center_of_mass())

    # random translations
    z_translate = slab_z_max - hdpe_copy.positions[:,2].min() + 3.0 + np.random.rand() * 2.0
    x_translate = (cell_x - chain_length) * 0.5 + (np.random.rand() - 0.5) * 0.2 * cell_x
    y_translate = cell_y * 0.5 + (np.random.rand() - 0.5) * 0.4 * cell_y
    hdpe_copy.translate([x_translate, y_translate, z_translate])

    combined = slab_copy + hdpe_copy
    combined.set_pbc([True, True, False])

    output_path = os.path.join(output_dir, f'config_{i:03d}.xyz')
    write(output_path, combined)
    



















