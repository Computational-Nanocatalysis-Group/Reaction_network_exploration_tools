#!/Users/user/opt/anaconda3/envs/qc/bin/python
# Don't forget to define the path to Python in your Anaconda environment

"""
Generating input for DIMER calculation
This script generates VASP input files (INCAR, POTCAR, KPOINTS, POSCAR) 
from an input structure (init.xyz by default). Use Vib2MODECAR.py for the generation of MODECAR.

Customize the parallelization parameters, NPAR and NCORE, to optimize performance for your machine!
See https://www.vasp.at/wiki/index.php/NPAR

You can now invoke the script with a file name argument:
    /path/to/this/script file_name.xyz
If no argument is provided, it defaults to 'init.xyz'.
"""

import sys
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase import Atoms
from typing import List
import os
from os import path

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Define simulation cell dimensions from geometry optimizations
cell = [20.8, 20.9, 21.0]  # e.g., Cu79
# cell = [21.4, 21.5, 21.6]  # e.g., Pd79

# Select metal
metal = 'Cu'
# metal = 'Pd'

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

def image_prep(atoms: Atoms, cell_par: List[float]) -> None:
    """
    Prepare the atomic image by setting the simulation cell, centering, 
    and enabling periodic boundary conditions.
    
    Parameters:
        atoms (Atoms): The atomic structure.
        cell_par (List[float]): List containing cell dimensions.
    """
    atoms.set_cell(cell_par)
    atoms.center()
    atoms.pbc = (True, True, True)


def sort_atoms(image: Atoms) -> Atoms:
    """
    VASP requires atoms to be sorted into element groups.
    
    This function sorts atoms in the image based on their element symbols.
    Atoms matching the 'metal' variable are given top priority.
    
    Parameters:
        image (Atoms): The atomic structure to be sorted.
        
    Returns:
        Atoms: A new Atoms object with metal atoms first, followed by the others sorted alphabetically.
    """
    sorted_indices = sorted(
        range(len(image)),
        key=lambda i: (0 if image[i].symbol == metal else 1, image[i].symbol)
    )
    return image[sorted_indices]

# -----------------------------------------------------------------------------
# VASP Calculator Setup
# -----------------------------------------------------------------------------

# Initialize the VASP calculator with parameters
calc_vasp = Vasp(
    prec='Normal',
    xc='revpbe',
    ivdw=12,
    encut=415.0,
    ediff=1.0E-06,
    algo='Normal',
    ismear=-1,
    sigma=0.03,
    nelm=400,
    lreal=False,
    iwavpr=11,
    ispin=2,
    isym=-1,
    lwave=False,
    lcharg=False,
    lasph=True,
    npar=2,  # Needs customization for your machine
    ncore=64,  # Needs customization for your machine
    icharg=2,
    ibrion=3,
    nsw=1000,
    lplane=True,
    potim=0,
    ediffg=-0.03,
    iopt=7,
    ichain=2
)

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

def main():
    # Use the file name provided as an argument; default to 'init.xyz' if none is given.
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'init.xyz'
    system = read(input_file)
    
    # Prepare the system: set cell, center, and periodic boundary conditions
    image_prep(system, cell)
    
    # Optionally, sort atoms by element (if needed)
    system = sort_atoms(system)
    
    # Set the VASP calculator and initialize
    system.calc = calc_vasp
    calc_vasp.initialize(system)
    
    # Uncomment the following lines to fix atoms of the selected metal if needed:
    # indices_to_fix = [atom.index for atom in system if atom.symbol == metal]
    # constraint = FixAtoms(indices=indices_to_fix)
    # system.set_constraint(constraint)
    
    # Write VASP input files: INCAR, POTCAR, KPOINTS, and POSCAR
    calc_vasp.write_incar(system)
    calc_vasp.write_potcar()
    calc_vasp.write_kpoints()
    write(images=system, filename='POSCAR', format="vasp")


if __name__ == '__main__':
    main()
