import os
import sys
import ase.io.lammpsdata as ld
import ase.io.extxyz as extxyz

phaseOffsets = {
    "amorph": 1,
    "divacancy": 39,
    "vacancy": 90,
    "interstitial": 60,
    "liquid": 123,
}

directory = os.getcwd()
data = "/../GAP_fitting/Training_Data/baseStructures.xyz"
inputFileName = directory + data
inputfile = open(inputFileName)

phase = sys.argv[1]
numberFiles = int(sys.argv[2])
start = int(sys.argv[3])
operation = sys.argv[4]
if operation.lower() == "lowanneal":
    operation = "lowAnneal"
elif operation.lower() == "highanneal":
    operation = "highAnneal"
elif operation.lower() == "medanneal":
    operation = "medAnneal"
elif operation.lower() == "heating":
    operation = "heating"
elif operation.lower() == "quenching":
    operation = "quenching"
elif not operation.lower() == "optimize":
    print("3rd argument must be one of: [optimize, lowAnneal, medAnneal, highAnneal]")
    print("Please try again with the proper argument")
    sys.exit()
iteration = sys.argv[5]

phaseStart = phaseOffsets[phase] + 1 + start
phaseEnd = phaseStart + numberFiles

for frame in range(phaseStart, phaseEnd):
    extxyzFrame = list(extxyz.read_extxyz(inputfile, frame))

    phase = extxyzFrame[0].info["config_type"]
    print(phase)

    out = (
        "/Si:H_LAMMPS_inputs/round"
        + str(iteration)
        + "/"
        + "GAP_"
        + operation
        + "_"
        + phase
        + "_"
        + str(frame - phaseOffsets[phase])
        + ".xyz"
    )
    outputFileName = directory + out
    outputfile = open(outputFileName, "w+")
    ld.write_lammps_data(outputfile, extxyzFrame)
