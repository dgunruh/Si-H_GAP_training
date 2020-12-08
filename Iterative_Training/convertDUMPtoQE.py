import ase.io.lammpsrun as lmprun
from collections import deque
import numpy as np
from ase.atoms import Atoms
import os
import ase.io.espresso as espresso
import ase.io.extxyz as extxyz
import shlex
import math
import sys

iteration = sys.argv[1]  # This tells the program which iteration round we are in
numFrames = int(sys.argv[2])

# specify any control parameters here
dt = "30.D0"
# calculation = "\"relax\""
calculation = ""
nstep = ""
pseudo_dir = "'./'"
outdir = "'./'"
pseudo_potential_Si = "Si.pbe-hgh.UPF"
pseudo_potential_H = "H.pbe-hgh.UPF"
occupation = '"smearing"'
smearing = '"marzari-vanderbilt"'
gauss = "0.05D0"
ecutwfc = "42.D0"

# convert the elements into their proper chemical symbols
def elm_conv(fld):
    return "14" if float(fld) == 2 else fld


def read_lammps_dump_text(fileobj, index=-1, **kwargs):
    """Process cleartext lammps dumpfiles

    :param fileobj: filestream providing the trajectory data
    :param index: integer or slice object (default: get the last timestep)
    :returns: list of Atoms objects
    :rtype: list
    """
    # Load all dumped timesteps into memory simultaneously
    lines = deque(fileobj.readlines())

    index_end = lmprun.get_max_index(index)

    n_atoms = 0
    images = []

    while len(lines) > n_atoms:
        line = lines.popleft()

        if "ITEM: TIMESTEP" in line:
            n_atoms = 0
            line = lines.popleft()
            # !TODO: pyflakes complains about this line -> do something
            # ntimestep = int(line.split()[0])  # NOQA

        if "ITEM: NUMBER OF ATOMS" in line:
            line = lines.popleft()
            n_atoms = int(line.split()[0])

        if "ITEM: BOX BOUNDS" in line:
            # save labels behind "ITEM: BOX BOUNDS" in triclinic case
            # (>=lammps-7Jul09)
            # !TODO: handle periodic boundary conditions in tilt_items
            tilt_items = line.split()[3:]
            celldatarows = [lines.popleft() for _ in range(3)]
            celldata = np.loadtxt(celldatarows)
            diagdisp = celldata[:, :2].reshape(6, 1).flatten()

            # determine cell tilt (triclinic case!)
            if len(celldata[0]) > 2:
                # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS"
                # to assign tilt (vector) elements ...
                offdiag = celldata[:, 2]
                # ... otherwise assume default order in 3rd column
                # (if the latter was present)
                if len(tilt_items) >= 3:
                    sort_index = [tilt_items.index(i) for i in ["xy", "xz", "yz"]]
                    offdiag = offdiag[sort_index]
            else:
                offdiag = (0.0,) * 3

            cell, celldisp = lmprun.construct_cell(diagdisp, offdiag)

            # Handle pbc conditions
            if len(tilt_items) > 3:
                pbc = ["pp" in d.lower() for d in tilt_items[3:]]
            else:
                pbc = (False,) * 3

        if "ITEM: ATOMS" in line:
            colnames = line.split()[2:]
            datarows = [lines.popleft() for _ in range(n_atoms)]
            data = np.loadtxt(datarows, converters={1: elm_conv})
            # print(data)
            out_atoms = lmprun.lammps_data_to_ase_atoms(
                data=data,
                colnames=colnames,
                cell=cell,
                celldisp=celldisp,
                atomsobj=Atoms,
                pbc=pbc,
                **kwargs
            )
            images.append(out_atoms)

        if len(images) > index_end >= 0:
            break

    return images[index]


def writeQE(input, output, phase, index=""):
    prefix = "'{}heatingIteration{}'".format(phase, index)

    # Now start to read in the LAMMPS file, and create the &SYSTEM namelist
    lines = open(input, "r").read().splitlines()
    line_iterator = 0
    file_iterator = 1
    nAtoms = 1
    infoDict = {}
    atomsXYZ = []
    nAtoms = 64
    for n, line in enumerate(lines):
        # get the number of atoms in the structure
        if line_iterator == 0:
            nAtoms = int(line)
        # Read in the header information, and put it into a dictionary
        elif line_iterator == 1:
            information = shlex.split(line)
            for j in information:
                fields = j.split("=")
                fieldname = fields[0]
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
            chemIndice = infoDict["Properties"]["species"][0]
            xyz_start = infoDict["Properties"]["pos"][0]
            xyz_end = infoDict["Properties"]["pos"][1]
            atomInfo = line.split()
            xyz = [float(s) for s in atomInfo[xyz_start:xyz_end]]
            xyz.append(atomInfo[chemIndice])
            atomsXYZ.append(xyz)
        line_iterator += 1

        if line_iterator > nAtoms + 1:
            print("File number: ", file_iterator)

            # We have reached the end of the structure, so now we need to do 2 things:
            # First, we need to create our map of bond lengths
            # Second, we need to add a specified number of H to our Si structure
            # Third, we need to create the DFT input file
            # Let's do the first task:
            latticeComponents = infoDict["Lattice"].split()
            latticeVectors = np.reshape(latticeComponents, (3, 3))

            # Get k-point numbers
            a1 = np.asarray([float(i) for i in latticeVectors[0]])
            a2 = np.asarray([float(i) for i in latticeVectors[1]])
            a3 = np.asarray([float(i) for i in latticeVectors[2]])

            a1Norm = np.linalg.norm(a1)
            a2Norm = np.linalg.norm(a2)
            a3Norm = np.linalg.norm(a3)

            k1 = math.ceil(2 * math.pi / (a1Norm * 0.2))
            k2 = math.ceil(2 * math.pi / (a2Norm * 0.2))
            k3 = math.ceil(2 * math.pi / (a3Norm * 0.2))

            # Now we need to create the DFT input file
            # First create the &CONTROL namelist
            out = open(output, "w")
            out.write("&CONTROL\n")

            if calculation:
                out.write("  calculation = " + calculation + ",\n")
            if dt:
                out.write("  dt = " + dt + ",\n")
            if prefix:
                out.write("  prefix=" + prefix + ",\n")
            if nstep:
                out.write("  nstep = " + nstep + ",\n")
            if pseudo_dir:
                out.write("  pseudo_dir = " + pseudo_dir + ",\n")
            if outdir:
                out.write("  outdir = " + outdir + "\n")
            out.write("  wf_collect = .true.,\n")
            out.write("  verbosity = 'high',\n")
            out.write("  tprnfor = .true.,\n")
            out.write("  tstress = .true.,\n")
            out.write("/\n")
            out.write("\n")

            # Now create the &SYSTEM namelist
            out.write("&SYSTEM\n")
            out.write("  ibrav = 0,")
            out.write(" nat = " + str(nAtoms) + ",")
            out.write(" ntyp = 2,\n")
            out.write("  ecutwfc\t= " + ecutwfc + ",\n")
            out.write("  degauss\t= " + gauss + ",\n")
            out.write("  occupations\t= " + occupation + ",\n")
            out.write("  smearing\t= " + smearing + ",\n")
            out.write("/\n")

            # Now create the &ELECTRONS namelist
            out.write("&ELECTRONS\n")
            out.write("  electron_maxstep = 200,\n")
            out.write("  conv_thr = 1.D-8,\n")
            out.write("  mixing_beta = 0.3D0,\n")
            out.write("  scf_must_converge=.false.,\n")
            out.write("/\n")

            # Now create the &IONS namelist
            out.write("&IONS\n")
            out.write("/\n")

            # Now create the ATOMIC_SPECIES card
            out.write("ATOMIC_SPECIES\n")
            out.write("Si  28.08  " + pseudo_potential_Si + "\n")
            out.write("H  1.00784  " + pseudo_potential_H + "\n")

            # Now create the CELL_PARAMETERS card
            out.write("CELL_PARAMETERS { angstrom }\n")
            for i in latticeVectors:
                out.write("   " + i[0] + "  " + i[1] + "  " + i[2] + "\n")

            out.write("K_POINTS {automatic}\n")
            out.write(str(k1) + " " + str(k2) + " " + str(k3) + "  0 0 0\n")

            out.write("ATOMIC_POSITIONS (angstrom)\n")
            for i in atomsXYZ:
                if i[3] == "Si":
                    out.write(
                        "Si "
                        + "{:.10f}".format(i[0])
                        + " "
                        + "{:.10f}".format(i[1])
                        + " "
                        + "{:.10f}".format(i[2])
                        + "\n"
                    )
                elif i[3] == "H":
                    out.write(
                        "H "
                        + "{:.10f}".format(i[0])
                        + " "
                        + "{:.10f}".format(i[1])
                        + " "
                        + "{:.10f}".format(i[2])
                        + "\n"
                    )


directory = os.getcwd()
lammpsDumpFolder = "/Si:H_LAMMPS_dumps/round" + str(iteration) + "/"
extxyzFolder = "/Si:H_LAMMPS_extxyzs/round" + str(iteration) + "/"
QEfolder = "/../Iterative_Training_DFT/Inputs/round" + str(iteration) + "/"

if not os.path.exists(directory + extxyzFolder):
    os.makedirs(directory + extxyzFolder)
if not os.path.exists(directory + QEfolder):
    os.makedirs(directory + QEfolder)

for inputfilename in os.listdir(directory + lammpsDumpFolder):
    inputfile = open(directory + lammpsDumpFolder + inputfilename)

    if numFrames == 1:
        # Read in the LAMMPS dump file
        images = read_lammps_dump_text(inputfile, index=-1)

        # Now get the output extxyz filename
        outputfile = directory + extxyzFolder + inputfilename[:-3] + "extxyz"

        begin = 4
        end = inputfilename.find("_", begin)
        operation = inputfilename[begin:end]

        begin = begin + len(operation) + 1
        end = inputfilename.find("_", begin)
        phase = inputfilename[begin:end]

        begin = begin + len(phase) + 5
        end = inputfilename.find(".xyz")
        number = inputfilename[begin:end]

        # for atoms in images:
        #    atoms.info = {"config_type":phase}

        # Output the extxyz file
        extxyz.write_extxyz(
            outputfile, images, columns=["symbols", "positions", "numbers"]
        )

        # Now make the QE input file
        qeName = phase + "_" + operation + number + ".in"
        qefile = directory + QEfolder + qeName
        writeQE(outputfile, qefile, phase)
    else:
        # Read in the LAMMPS dump file
        images = read_lammps_dump_text(inputfile, index=slice(1, numFrames + 1, 1))
        for i, image in enumerate(images):
            # Now get the output extxyz filename
            outputfile = (
                directory
                + extxyzFolder
                + inputfilename[:-4]
                + "_"
                + str(i + 1)
                + ".extxyz"
            )

            begin = 4
            end = inputfilename.find("_", begin)
            operation = inputfilename[begin:end]

            begin = begin + len(operation) + 1
            end = inputfilename.find("_", begin)
            phase = inputfilename[begin:end]

            begin = begin + len(phase) + 11
            end = inputfilename.find(".xyz")
            number = inputfilename[begin:end]

            # for atoms in images:
            #    atoms.info = {"config_type":phase}

            # Output the extxyz file
            extxyz.write_extxyz(
                outputfile, image, columns=["symbols", "positions", "numbers"]
            )

            # Now make the QE input file
            qeName = phase + "_" + operation + number + "_" + str(i + 1) + ".in"
            qefile = directory + QEfolder + qeName
            writeQE(outputfile, qefile, phase, str(i + 1))
