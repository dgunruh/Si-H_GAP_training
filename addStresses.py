#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import shlex
import math
import pipes

# Get the input filename
inputFileName = sys.argv[1]
outputFileName = inputFileName[:-4] + "_addedStresses.xyz"

# First get the file containing all the structures
directory = os.getcwd()
file_extension = "/GAP_fitting/Training_Data/"
inputfile = directory + file_extension + inputFileName

# Create the output file
outputfile = directory + file_extension + outputFileName

# Now create the &CONTROL namelist
out = open(outputfile, "w")


# Now read in the extended xyz file
lines = open(inputfile, "r").read().splitlines()
line_iterator = 0
file_iterator = 1
nAtoms = 1
infoDict = {}
atomsExtended = []
atomsForces = []
e_sigma = 0.001
f_sigma = 0.1
s_sigma = 0.05
v_sigma = 0
information = []
atomLines = []
for n, line in enumerate(lines):
    # get the number of atoms in the structure
    if line_iterator == 0:
        nAtoms = int(line)
    # Read in the header information, and put it into a dictionary
    elif line_iterator == 1:
        quoteInformation = []
        # quotedinformation = shlex.split(line, posix=False)
        information = shlex.split(line)
        for j in information:
            fields = j.split("=")
            fieldname = fields[0]
            splitfield = fields[1].split(" ")
            print(fields)
            if len(splitfield) >= 2:
                quoteField = ""
                for n, k in enumerate(splitfield):
                    quote = k + " "
                    if n == 0:
                        quote = '"' + k + " "
                    if n == len(splitfield) - 1:
                        quote = k + '"'
                    quoteField += quote
                quoteInformation.append(str(fields[0] + "=" + quoteField))
            else:
                quoteInformation.append(fields[0] + "=" + fields[1])

            if fieldname == "Properties":
                properties = fields[1].split(":")
                fielddata = {}
                start = 0
                for n, l in enumerate(properties):
                    if n % 3 == 0:
                        item = l
                    elif n % 3 == 2:
                        end = start + int(l)
                        fielddata[item] = [start, end]
                        start = end
                infoDict[fieldname] = fielddata
            else:
                infoDict[fieldname] = fields[1]
    # Read in the extended xyz atom data
    elif line_iterator <= nAtoms + 1:
        atomLines.append(line)
        atomInfo = line.split()
        force_start = infoDict["Properties"]["dft_force"][0]
        force_end = infoDict["Properties"]["dft_force"][1]
        for s in atomInfo[force_start:force_end]:
            atomsForces.append(abs(float(s)))
    line_iterator += 1

    if line_iterator > nAtoms + 1:
        print("File number: ", file_iterator)

        # We have reached the end of the structure, so now we need to do 2 things:
        stresses = []
        stressString = infoDict["dft_virial"]
        # stressSplit = stressString.split(" ")
        for l in stressString.split(" "):
            stresses.append(abs(float(l)))
        mean_stress = np.mean(np.array(stresses))
        s_sigma = abs(0.1 * mean_stress)

        mean_force = np.mean(np.array(atomsForces))
        if mean_force >= 0.1:
            f_sigma = abs(0.1 * mean_force)
        else:
            f_sigma = 0.01

        sigma_string = (
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

        quoteInformation.insert(1, sigma_string)
        newline = " ".join(quoteInformation)

        out.write(str(nAtoms) + "\n")
        out.write(newline + "\n")
        for aline in atomLines:
            out.write(aline + "\n")

        # reset the iterator
        line_iterator = 0
        file_iterator += 1
        atomsXYZ = []
        atomsForces = []
        atomLines = []
