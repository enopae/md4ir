"""
Extract variables from xyz data
"""

__all__ = [ 'extract' ]

import numpy as np
import csv
from .fileHandling import *

def getCentroid(coords):
    centroid_coords = []
    for i in range(len(coords[0])):
        xyz = np.zeros(3)
        for j in range(3):
            xyz[j] = coords[0][i][j] + (coords[1][i][j] - coords[0][i][j])/2
        centroid_coords.append(xyz)
    return centroid_coords



class _Extract:
    """Extracts the coordinates into a coords object
    """

    def __init__(self,atoms,file):
        """Initialize
        
        Arguments:
            atoms {list} -- list of atom indices
            file {string} -- xyz file name
        """
        self.atoms = atoms          # list of atoms
        self.data = read_file(file) # data from the file
        self.natoms = int(len(self.atoms)) # number of atoms
        # Coordinates for each atom in the xyz file
        self.coords = [[] for i in range(self.natoms)]
        for i in range(len(self.coords)):
            self.coords[i] = self.coord2float(int(self.atoms[i]))

    def coord2float(self,idx):
        """Converts the atomic coordinates of atom number idx to floats from xyz data.
                
        Arguments:
            idx {int} -- atom index
        """
        length = int(self.data[0][0])
        length += 2 # distance from one atom coordinate to the next
        coords = []
        idx = idx + 1 # atom index in python (including header)
        for i, line in enumerate(self.data):
            if i == idx:
                xyz = line[1:]
                xyz = [float(c) for c in xyz]
                coords.append(xyz)
                idx += length
        return coords

class _Calc:
    
    def __init__(self,coords):
        """Initialize
        
        Arguments:
            coords {list} -- list of lists for the coordinates of each atom
        """
        self.coords = coords
        self.natoms = int(len(self.coords)) # number of atoms
        self.calculators = { # calculators for different number of atoms
                            2: self.calc_bond,
                            3: self.calc_angle,
                            4: self.calc_dihedral}
        self.datalength = len(self.coords[0])
        # Calculate the values
        self.values = []
        self.calculators[self.natoms]()

    def calc_bond(self):
        """Calculates the bond lengths
        """
        for i in range(self.datalength):
            d = 0
            print(self.coords[0][i],self.coords[1][i])
            for j in range(3):
                d += (self.coords[0][i][j] - self.coords[1][i][j])**2
            d = np.sqrt(d)
            self.values.append(d)

    def calc_angle(self):
        """Calculates the bond angles
        """
        for i in range(self.datalength):
            # Define bond vectors
            v1 = np.zeros(3)
            v2 = np.zeros(3)
            for j in range(3):
                v1[j] = self.coords[0][i][j] - self.coords[1][i][j]
                v2[j] = self.coords[2][i][j] - self.coords[1][i][j]
            cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            self.values.append(np.arccos(cos)*180/np.pi) # append and convert to degrees

    def calc_dihedral(self):
        """Calculates the dihedral angles
        """
        # Algorithm is based on:
        # https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
        for i in range(self.datalength):
            # Define bond vectors
            v1 = np.zeros(3)
            v2 = np.zeros(3)
            v3 = np.zeros(3)
            for j in range(3):
                v1[j] = self.coords[0][i][j] - self.coords[1][i][j]
                v2[j] = self.coords[2][i][j] - self.coords[1][i][j]
                v3[j] = self.coords[3][i][j] - self.coords[2][i][j]
            v2 /= np.linalg.norm(v2)
            tmp1 = v1 - np.dot(v1, v2)*v2
            tmp2 = v3 - np.dot(v2, v3)*v2
            x = np.dot(tmp1, tmp2)
            y = np.dot(np.cross(v2, tmp1), tmp2)
            dihedral_angle = np.degrees(np.arctan2(y,x))
            self.values.append(dihedral_angle)


"""
To be rechecked...
    def loop_bonds(bond_list, sum = True, scale = 1.0, timestep = 0.2, t_start = 0, t_end = 20, wn_start = 0, wn_end = 4000, trajectory_file = 'xtb.trj'):
        signal = 0
        for atom1,atom2 in bond_list:
            bonds = get_bonds(atom1, atom2, trajectory_file)
            bond = Calc(bonds, scale, timestep, t_start, t_end, wn_start, wn_end)
            if sum == True:
                bond.gen_spectrum(norm = False)
                signal = signal + bond.signal
        signal = signal / np.amax(signal)
        return bond.wavenumbers, signal

    def loop_angles(angle_list, sum = True, scale = 1.0, timestep = 0.2, t_start = 0, t_end = 20, wn_start = 0, wn_end = 4000, trajectory_file = 'xtb.trj'):
        signal = 0
        for atom1,atom2,atom3 in angle_list:
            angles = get_angles(atom1, atom2, atom3, trajectory_file)
            angle = Calc(angles, scale, timestep, t_start, t_end, wn_start, wn_end)
            if sum == True:
                angle.gen_spectrum(norm = False)
                signal = signal + angle.signal
        signal = signal / np.amax(signal)
        return angle.wavenumbers, signal
"""
# Execute
        
def extract(args):
   
    # Set up variables
    name = args.name
    file = args.file
    atoms = args.atoms.split(',')

    # Run routines
    if args.centroid:
        centroid = args.centroid.split(',')
        extractor = _Extract(atoms, file)
        centroid_extractor = _Extract(centroid, file)
        centroid_coords = getCentroid(centroid_extractor.coords)
        coords = [extractor.coords[0], centroid_coords]
        calculator = _Calc(coords)
    else:
        extractor = _Extract(atoms, file)
        calculator = _Calc(extractor.coords)
    
    data = calculator.values

    # Write the data into a file
    with open(name+'_'+args.atoms.replace(',', '-')+'.txt','w',encoding='utf8') as f:
        for i in np.arange(0,len(data)):
                f.write("%f\n" % (data[i]))
        f.close()
    return