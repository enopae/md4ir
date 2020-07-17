#!/bin/env python3
"""
md4ir
MD analysis for IR spectra calculations
Version: 0.1
Written by: Eno Paenurk
"""

import md4ir # local modules
import argparse as ap

def main():
    # Parse arguments
    parser = ap.ArgumentParser()
    
    # Set up parent parser for common arguments
    parent_parser = ap.ArgumentParser(add_help=False)
    parent_parser.add_argument("-n", "--name", help="Basename for the output file.", type=str, default='dummy')

    # Set up subparsers for modules
    subparsers = parser.add_subparsers(help = 'Available methods.')

    # Parser for calculateSpectra
    spec_parser = subparsers.add_parser('spec', help = 'Calculate the spectrum.', formatter_class=ap.ArgumentDefaultsHelpFormatter, parents=[parent_parser])
    spec_parser.set_defaults(func=md4ir.calcSpec)
    spec_parser.add_argument("-f", "--file", 
        help="The data file (dipoles, bond lengths, etc). Loop over files: use -f for each file. Average: use -f on a comma-separated list of files.",
        action='append', type=str)
    spec_parser.add_argument("-s", "--scale", help="Scaling factor for the wavenumbers of the computed spectrum.", type=float, default=1.0)
    spec_parser.add_argument("--timestep", help="Simulation timestep in fs.", type=float, default=0.5)
    spec_parser.add_argument("--time_start", help="Start evaluating the simulation at that time (in ps).", type=int, default=0)
    spec_parser.add_argument("--time_end", help="End evaluating the simulation at that time (in ps).", type=int, default=20)
    spec_parser.add_argument("--wn_start", help="Start evaluating the simulation at that wavenumber (in cm-1).", type=int, default=0)
    spec_parser.add_argument("--wn_end", help="End evaluating the simulation at that wavenumber (in cm-1).", type=int, default=4000)
    spec_parser.add_argument("--method", help="Spectrum calculation method.", choices = ['numpy'], type=str, default = 'numpy')

    # Parser for extractVariables
    extract_parser = subparsers.add_parser('extract', help='Extract variables from an xyz trajectory.', parents=[parent_parser])
    extract_parser.set_defaults(func=md4ir.extract)
    extract_parser.add_argument("-f", "--file", help="Trajectory file.", type=str, default='traj.xyz')
    extract_parser.add_argument("-a", "--atoms", help="String of atoms defining the variable, comma separated (2=bond, 3=angle, 4=dihedral).")
    extract_parser.add_argument("-c", "--centroid", help="String of two atoms for the centroid, comma separated.")

    # Parser for calculateDipoles
    dipole_parser = subparsers.add_parser('dipole', help = 'Calculate dipoles from charges and an xyz trajectory.', parents=[parent_parser])
    dipole_parser.set_defaults(func=md4ir.calcDipole)
    dipole_parser.add_argument("-f", "--file", help="Trajectory file.", type=str, default='traj.xyz')
    dipole_parser.add_argument("-c", "--charges", help="File for list of charges.", type=str, default='charges')

    # Collect the arguments
    args = parser.parse_args()

    # Execute the function
    args.func(args)

if __name__ == "__main__":
    main()
