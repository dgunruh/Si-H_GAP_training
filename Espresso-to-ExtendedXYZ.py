#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import shlex
import math

# Define constants
bohr_to_angstrom = 0.529177
ry_to_ev = 13.605698
config_type = sys.argv[4]
iterationRound = sys.argv[1]

# Define default sigma values
e_sigma = 0.001
f_sigma = 0.1
s_sigma = 0.05
v_sigma = 0


# First get the QE output file (in the verbose format)
directory = os.getcwd()
inputFileName = sys.argv[2]
# inputfile = directory + "/" +  "Si:H_DFT_structure_outputs" + "/Hydrogenated_Si_Structures-round2/" + inputFileName
inputfile = (
    directory
    + "/"
    + "Iterative_Training_DFT/Outputs/round"
    + iterationRound
    + "/"
    + inputFileName
)

lines = open(inputfile, "r").read().splitlines()

alat = 0.0
volume = 0.0
nAtoms = 0
axes = [[], [], []]
energy = 0.0
stress = []
atoms = []
forces = []

crystalAxesEnd = 0
coordsStart = 0
forceStart = 0
stressStart = 0
getCoords = False
getForces = False
getEnergy = False
getStress = False

forceComponents = []
stressComponents = []

for n, line in enumerate(lines):
    newline = line.split()
    # print(newline)
    # print(len(newline))
    if len(newline) > 2:
        if newline[0] == "unit-cell":
            volume = float(newline[3]) * (bohr_to_angstrom ** 3)
        if newline[2] == "(alat)":
            alat = float(newline[4]) * bohr_to_angstrom
        if newline[2] == "atoms/cell":
            nAtoms = int(newline[4])
        if newline[0] == "crystal":
            crystalAxesEnd = n + 4
        if n < crystalAxesEnd and n >= crystalAxesEnd - 3:
            ax0 = float(newline[3])
            ax1 = float(newline[4])
            ax2 = float(newline[5])
            axvec = [str(ax0 * alat), str(ax1 * alat), str(ax2 * alat)]
            axes[n - (crystalAxesEnd - 3)] = axvec

        if newline[1] == "Fermi":
            getEnergy = True

        if (
            getEnergy
            and newline[0] == "!"
            and newline[1] == "total"
            and newline[2] == "energy"
        ):
            energy = float(newline[4]) * ry_to_ev

        if newline[0] == "total" and newline[1] == "stress":
            getStress = True
            stressStart = n

    if len(newline) > 4:
        if newline[4] == "(cartesian":
            forceStart = n + 1
            getForces = True
    if line == "   Cartesian axes":
        coordsStart = n + 2
        getCoords = True
    if n > coordsStart and n < coordsStart + nAtoms + 1 and getCoords:
        atoms.append(
            [
                newline[1],
                str(float(newline[6]) * alat),
                str(float(newline[7]) * alat),
                str(float(newline[8]) * alat),
            ]
        )
        if newline[1] == "H":
            atoms[n - coordsStart - 1].append("1")
        else:
            atoms[n - coordsStart - 1].append("14")
    if n > forceStart and n < forceStart + nAtoms + 1 and getForces:
        atoms[n - forceStart - 1].extend(
            [
                str(float(newline[6]) * ry_to_ev / bohr_to_angstrom),
                str(float(newline[7]) * ry_to_ev / bohr_to_angstrom),
                str(float(newline[8]) * ry_to_ev / bohr_to_angstrom),
            ]
        )
        force1 = float(newline[6]) * ry_to_ev / bohr_to_angstrom
        force2 = float(newline[7]) * ry_to_ev / bohr_to_angstrom
        force3 = float(newline[8]) * ry_to_ev / bohr_to_angstrom
        if math.sqrt(force1 ** 2 + force2 ** 2 + force3 ** 2) > 15:
            print(
                "Error: force on atom too large, greater than 15 eV/angstrom. Filename is: ",
                inputFileName,
            )
        forceComponents.append(force1)
        forceComponents.append(force2)
        forceComponents.append(force3)

    if n > stressStart and getStress and n <= stressStart + 3:
        stress.append(
            [
                str(float(newline[0]) * ry_to_ev / (bohr_to_angstrom ** 3) * volume),
                str(float(newline[1]) * volume * ry_to_ev / (bohr_to_angstrom ** 3)),
                str(float(newline[2]) * volume * ry_to_ev / (bohr_to_angstrom ** 3)),
            ]
        )

        stress1 = float(newline[0]) * ry_to_ev / (bohr_to_angstrom ** 3) * volume
        stress2 = float(newline[1]) * volume * ry_to_ev / (bohr_to_angstrom ** 3)
        stress3 = float(newline[2]) * volume * ry_to_ev / (bohr_to_angstrom ** 3)
        stressComponents.append(stress1)
        stressComponents.append(stress2)
        stressComponents.append(stress3)

if config_type == "amorph":
    e_sigma = 0.01
elif config_type == "liquid":
    e_sigma = 0.003

forceArray = np.array(forceComponents)
fmean = np.mean(forceArray)
if fmean >= 0.1:
    f_sigma = 0.1 * np.mean(forceArray)
else:
    fmean = 0.01

stressArray = np.array(stressComponents)
s_sigma = 0.1 * np.mean(stressArray)


# Now output the extended xyz file
outputFileName = sys.argv[3]
outputFolder = (
    directory
    + "/"
    + "Si:H_Extended_XYZ_structures/Iterative_Training_round"
    + iterationRound
    + "/"
)

if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
outputfile = outputFolder + outputFileName
out = open(outputfile, "w")
out.write(str(nAtoms) + "\n")

# Create the properties string
config = "config_type=" + config_type
sigma = (
    'sigma="'
    + str(e_sigma)
    + " "
    + str(f_sigma)
    + " "
    + str(s_sigma)
    + " "
    + str(v_sigma)
    + '"'
)
de = "dft_energy=" + str(energy)
virial = (
    "dft_virial="
    + '"'
    + " ".join(stress[0])
    + " "
    + " ".join(stress[1])
    + " "
    + " ".join(stress[2])
    + '"'
)

pbc = 'pbc="T T T"'
lattice = (
    "Lattice="
    + '"'
    + " ".join(axes[0])
    + " "
    + " ".join(axes[1])
    + " "
    + " ".join(axes[2])
    + '"'
)
properties = "Properties=species:S:1:pos:R:3:Z:I:1:dft_force:R:3"

descriptor = (
    config
    + " "
    + sigma
    + " "
    + de
    + " "
    + virial
    + " "
    + pbc
    + " "
    + lattice
    + " "
    + properties
    + "\n"
)
out.write(descriptor)

for i in atoms:
    out.write("\t".join(i) + "\n")
