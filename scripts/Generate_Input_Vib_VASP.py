#!/Users/mpol/opt/anaconda3/envs/qc/bin/python
# Don't forget to define the path to Python in your Anaconda environment
"""
Generate VASP input files for vibrational frequency calculations.

This script reads an input structure from CONTCAR and generates VASP input files
(INCAR, POTCAR, KPOINTS, POSCAR).

Customize the parallelization parameters, NPAR and NCORE, to optimize performance for your machine!
See https://www.vasp.at/wiki/index.php/NPAR

Usage (no file selected by default):
    /path/to/this/script

"""

import os
from os import path
from typing import List

from ase.io import read, write
from ase.io.vasp import read_vasp
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase import Atoms

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Select metal for potential constraints
metal = 'Cu'
# metal = 'Pd'

# For DIMER calculations, may use 'CENTCAR'
input_file = 'CONTCAR'

# -----------------------------------------------------------------------------
# VASP Calculator Setup
# -----------------------------------------------------------------------------

calc_vasp = Vasp(
    prec='Normal',
    xc='revpbe',
    ivdw=12,
    encut=415.0,
    ediff=1.0E-07,
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
    npar=2, # Needs customization for your machine
    ncore=64,  # Needs customization for your machine
    icharg=1,
    ibrion=5,
    maxmix=60,
    nfree=2,
    lplane=True
)

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

def main():
    # Read the input structure from the specified VASP file (CONTCAR)
    system = read_vasp(input_file)
    
    # Set the VASP calculator for the system and initialize it
    system.calc = calc_vasp
    calc_vasp.initialize(system)
    
    # Apply constraints: fix all atoms of the selected metal
    indices_to_fix = [atom.index for atom in system if atom.symbol == metal]
    system.set_constraint(FixAtoms(indices=indices_to_fix))
    
    # Generate VASP input files: INCAR, POTCAR, KPOINTS, and POSCAR
    calc_vasp.write_incar(system)
    calc_vasp.write_potcar()
    calc_vasp.write_kpoints()
    write('POSCAR', system, format="vasp")


if __name__ == '__main__':
    main()
