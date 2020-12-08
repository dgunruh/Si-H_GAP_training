#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import shlex
import math

#Define the two atoms
atom1 = "Si"
atom2 = "Si"

#specify any control parameters here
prefix = "'{}{}'".format(atom1, atom2)
dt = "30.D0"
#calculation = "\"relax\""
calculation = ""
nstep = ""
pseudo_dir = "'./'"
outdir = "'./'"
pseudo_potential_Si = "Si.pbe-hgh.UPF"
pseudo_potential_H = "H.pbe-hgh.UPF"
ecutwfc = "42.D0"
basisLength = 20.0
kpoints = math.ceil(2*math.pi/(basisLength*0.2))
#kpoints = 8

#Choose 2.7 angstroms as our Si-Si bond length
SiSibondLength = 2.5
#Choose 1.58 angstroms as our Si-H bond length
SiHbondLength = 1.58
#Choose 74 picometers as our H-H bond length
HHbondLength = .74

bond = 0.0
if atom1 == "Si" and atom2 == "Si":
	bond = SiSibondLength
elif atom1 == "Si" and atom2 == "H":
	bond = SiHbondLength
else:
	bond = HHbondLength

nPoints = 25
coarsePoints = 3
finerPoints = 2
fineClosePoints = 10
finestPoints = nPoints - coarsePoints - finerPoints - fineClosePoints - 1
stepSize = bond/100

coarseStart = bond - coarsePoints*stepSize*5
fineStart = coarseStart - finerPoints*stepSize*3
fineEnd = coarseStart - stepSize*3
finestStart = fineStart - finestPoints*stepSize*2
finestEnd = fineStart - stepSize
end = bond
fineStartClose = finestStart - (fineClosePoints)*stepSize*2.5
fineEndClose = finestStart - stepSize*2.5

finerLengthsClose = np.linspace(fineStartClose, fineEndClose, fineClosePoints)
finestLengths = np.linspace(finestStart, finestEnd, finestPoints)
finerLengths = np.linspace(fineStart, fineEnd, finerPoints)
coarseLengths = np.linspace(coarseStart, end, coarsePoints + 1)
bondLengths = np.concatenate((finestLengths, finerLengths, coarseLengths, finerLengthsClose))
print(bondLengths)

if bond < 1:
	roundedLengths = np.array([round(i,3) for i in bondLengths])
else:
	roundedLengths = np.array([round(i,2) for i in bondLengths])


if atom1 == atom2:
	nTypes = 1
else:
	nTypes = 2

#################################
# Now we actually run the code! #
#################################

#Get the output filefolder
directory = os.getcwd()
outputFileFolder = "/" + atom1 + "-" + atom2 + " DFT input files/"

#Write out a file containing the two atoms at each distance
for n, l in enumerate(roundedLengths):
	#Create and open the output file
	outputfile = directory + outputFileFolder + atom1 + atom2 + str(n + 1) + ".in"
	out = open(outputfile,'w')

	#Now create the &CONTROL namelist
	out.write("&CONTROL\n")
	if calculation:
		out.write("  calculation = "+calculation+ ",\n")
	if dt:
		out.write("  dt = " + dt + ",\n")
	if prefix:
		out.write("  prefix=" + prefix + ",\n")
	if nstep:
		out.write("  nstep = "+nstep+",\n")
	if pseudo_dir:
		out.write("  pseudo_dir = "+pseudo_dir+ ",\n")
	if outdir:
		out.write("  outdir = " + outdir+ "\n")
	out.write("  wf_collect = .true.,\n")
	out.write("  verbosity = 'high',\n")
	out.write("  tprnfor = .true.,\n")
	out.write("  tstress = .true.,\n")
	out.write("/\n")
	out.write("\n")

	#Now create the &SYSTEM namelist
	out.write("&SYSTEM\n")

	#Make a big enough unit cell that the atoms won't interact with themselves via pbc
	out.write("  ibrav = 1,")
	out.write(" A = " + str(basisLength) + ",")
	out.write(" nat = 2,")
	out.write(" ntyp = " + str(nTypes) + ",\n")
	out.write("  ecutwfc\t= " + ecutwfc + ",\n")
	out.write("  degauss\t= 0.05D0,\n")
	out.write("  occupations\t= \"smearing\",\n")
	out.write("  smearing\t= \"marzari-vanderbilt\",\n")
	out.write("/\n")
	out.write("\n")

	#Now create the &ELECTRONS namelist
	out.write("&ELECTRONS\n")
	out.write("  electron_maxstep = 200,\n")
	out.write("  conv_thr = 1.D-8,\n")
	out.write("  mixing_beta = 0.3D0,\n")
	#out.write("  scf_must_converge=.false.,\n")
	out.write("/\n")
	out.write("\n")

	#Now create the &IONS namelist
	#out.write("&IONS\n")
	#out.write("/\n")
	#out.write("\n")

	#Now create the ATOMIC_SPECIES card
	out.write("ATOMIC_SPECIES\n")
	if atom1 == "Si":
		out.write("Si  28.08  " + pseudo_potential_Si + "\n")
	if atom2 == "H":
		out.write("H  1.00784  " + pseudo_potential_H + "\n")
	out.write("\n")

	#Now create the CELL_PARAMETERS card
	#out.write("CELL_PARAMETERS { angstrom }\n")
	#for i in latticeVectors:
#		out.write("   " + i[0] + "  " + i[1] + "  " + i[2] + "\n")

	#Now output the K_POINTS
	out.write("K_POINTS {automatic}\n")
	out.write(str(kpoints) + " " + str(kpoints) + " " + str(kpoints) + "  0 0 0\n")
	out.write("\n")

	#Now output the atomic locations
	out.write("ATOMIC_POSITIONS (angstrom)\n")
	#put atom 1 at the origin
	out.write(atom1 + " {:.10f}".format(0.0) + " " + "{:.10f}".format(0.0) + " " + "{:.10f}".format(0.0) + "\n")
	out.write(atom2 + " {:.10f}".format(l) + " " + "{:.10f}".format(0.0) + " " + "{:.10f}".format(0.0) + "\n")
