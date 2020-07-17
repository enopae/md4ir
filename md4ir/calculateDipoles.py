"""
Calculate dipoles from charges and an xyz trajectory
"""

__all__ = [ 'calcDipole' ]

import numpy as np
from . import fileHandling

def calcDipole(name, file, charges):
    # Set up variables
    #name = args.name
    #datafile = args.file
    #charges = args.charges

    # Run routines
    
    # Calculate dipoles
    charge_list = get_charges(charges)
    atoms = get_atoms(datafile)
    dipoles = calc_dipole(charge_list, atoms)

    # Write the dipoles into a file
    with open(name+'_'+'dipoles.txt','w',encoding='utf8') as f:
        for i in np.arange(0,len(dipoles)):
                f.write("%f %f %f\n" % (dipoles[i][0],dipoles[i][1],dipoles[i][2]))
        f.close()
    return


# Functions

def get_charges(chargefile):
    lines = fileHandling.read_file(chargefile)
    charges = []
    for line in lines:
        charges.append(float(line[0]))
    return charges

def get_atoms(xyzfile):
    lines = fileHandling.read_file(xyzfile)
    n_atoms = int(lines[0][0])
    n_lines = n_atoms + 2
    n_xyzs = int(len(lines) / n_lines)
    data = [[[] for j in range(n_atoms)] for i in range(n_xyzs)]
    i=2
    for xyz in range(n_xyzs):
        for atom in range(0,n_atoms):
            data[xyz][atom] = lines[i+atom][1:]
        i += n_lines
    return data

def calc_dipole(charges, xyzs):
    dipoles = []
    for xyz in xyzs:
        charge_coord = [0,0,0]
        for atom in range(0,len(charges)):
            for i in range(0,3):
               charge_coord[i] += (charges[atom]*float(xyz[atom][i]))
        dipoles.append(charge_coord)
    return dipoles
        