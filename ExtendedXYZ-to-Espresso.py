#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import shlex
import math

#Define the phase of the pure Si
phase = sys.argv[1][:-4]
#phase = "diamondInterstitial"

#specify any control parameters here
prefix = "'{}SiH'".format(phase)
dt = "30.D0"
#calculation = "\"relax\""
calculation = ""
nstep = ""
pseudo_dir = "'./'"
outdir = "'./'"
pseudo_potential_Si = "Si.pbe-hgh.UPF"
pseudo_potential_H = "H.pbe-hgh.UPF"
occupation = "\"smearing\""
smearing = "\"marzari-vanderbilt\""
gauss = "0.05D0"
ecutwfc = "42.D0"
natoms_total = 64

#Choose 2.7 angstroms as our Si-Si bond length
SiSibondLength = 2.7
#Choose 1.58 angstroms as our Si-H bond length 
SiHbondLength0 = 1.58

SiRadius = 1.11
HRadius = 0.38

Sigma_bl = .1 #Angstroms
H2angle = 92 #degrees
percentInterstitial = 0.50
percentDihydride = .3
siUnitCellLength = 5.43

#The atomic concentration of H we are adding to the pure-Si structures
# (defined relative to the number of Si atoms)
percent_H = .12
# if phase == "liquid" or phase == "amorphous":
# 	percent_H = .12
# else:
# 	percent_H = .1

def getDistanceMap(atomXYZ, latticeInfo,bondLength):
	a1 = np.asarray([float(i) for i in latticeInfo[0]])
	a2 = np.asarray([float(i) for i in latticeInfo[1]])
	a3 = np.asarray([float(i) for i in latticeInfo[2]])

	#determine frac_pos for each atom if not already given (position in terms of a1,a2,a3)
	A = np.dstack((a1,a2,a3))
	A_inverse = np.linalg.inv(A)[0]
	frac_lattice_positions = []
	for i in atomXYZ:
		pos = np.asarray(i)
		frac_pos = np.dot(A_inverse,pos)
		frac_lattice_positions.append(frac_pos)

	#from frac_lattice_positions, we can determine periodic boundary conditions by comparing 
	#fractional values of lattice vectors
	bondedNeighbors = {}
	for n, i in enumerate(frac_lattice_positions):
		bondedNeighbors[n] = []
		#for j in np.concatenate(frac_lattice_positions[:n], frac_lattice_positions[n+1:]):
		for j in frac_lattice_positions[:n] + frac_lattice_positions[n+1:]:
			lattice_diff = np.zeros(3)
			for k in range(3):
				raw = j[k]-i[k]
				diff = raw
				#print(raw)
				if abs(raw) > 0.5:
					#use pbc
					if raw < 0:
						diff = 1 + raw
					else:
						diff = raw - 1
				lattice_diff[k] = diff
			#print(lattice_diff)
			delta = lattice_diff[0]*a1 + lattice_diff[1]*a2 + lattice_diff[2]*a3
			#print(delta)
			distance = np.linalg.norm(delta)
			if distance < bondLength:
				bondedNeighbors[n].append(delta)

	#break

	return bondedNeighbors

def getH_XYZ(atomsXYZ, bondMap, num_H):
	#Get all atoms which have dangling, floating and strained bonds
	dihydrideSites = 0
	dihydrideAtoms = []
	danglingBonds = 0
	danglingAtoms = []
	floatingBonds = 0
	floatingAtoms = []
	strainedBonds = 0
	strainedAtoms = []
	interstitialSites = 0
	interstitialAtoms = []
	bondangles = []
	if phase != "liquid":
		for n, bondVectors in bondMap.items():
			atom = atomsXYZ[n]
			#Make the bond length of Si-H have a Gaussian width (to increase diversity of atomic neighborhoods)
			SiHbondLength = np.random.normal(SiHbondLength0,Sigma_bl)
			
			#Get the number of neighbors to the Si atom
			numNeighbors = len(bondVectors)
			
			#SiH2 complex
			if numNeighbors == 2:
				dihydrideSites += 1

				Si_1 = bondVectors[0]
				Si_2 = bondVectors[1]
				interBondAngle = np.arccos(np.dot(Si_1, Si_2)/(np.linalg.norm(Si_1)*np.linalg.norm(Si_2)))
				theta = -((360 - interBondAngle*180/np.pi)/2 - H2angle/2)*np.pi/180

				z = np.cross(Si_1, Si_2)
				y = np.cross(Si_1,z)/np.linalg.norm(z)
				rot1 = np.cos(theta)*Si_1 + np.sin(theta)*y
				rot1 = rot1*SiHbondLength/np.linalg.norm(rot1)
				rot2 = np.cos(-theta+interBondAngle)*Si_1 + np.sin(-theta + interBondAngle)*y
				rot2 = rot2*SiHbondLength/np.linalg.norm(rot2)

				dihydrideAtoms.append([atom+rot1,atom + rot2])

			#dangling bond
			elif numNeighbors == 3:
				danglingBonds+=1

				#Get the dangling bond vector
				a1 = bondVectors[0]
				a2 = bondVectors[1]
				a3 = bondVectors[2]
				d = np.cross(a1,a2) + np.cross(a2,a3) + np.cross(a3,a1)
				d = SiHbondLength*d/np.linalg.norm(d)

				danglingAtoms.append([atom+d])

			#floating bond(s)
			elif numNeighbors >= 5:
				floatingBonds += 1

				added = False
				for j in bondVectors:
					if not added:
						f = -SiHbondLength*j/np.linalg.norm(j)
						
						dontAdd = False
						for i in bondVectors:
							separation = np.linalg.norm((atom + i) - (atom + f))
							if separation < (SiRadius + HRadius):
								dontAdd = True
						if not dontAdd:
							floatingAtoms.append([atom + f])
							added = True
				
				iterations = 0
				while not added:
					#Here there is no good spot to add H, except at random next to the central Si
					f0 = np.array([np.random.randint(-1,2),np.random.randint(-1,2),np.random.randint(-1,2)])
					fhat = f0/np.linalg.norm(f0)
					f = SiHbondLength*fhat
					dontAdd = False
					for i in bondVectors:
						separation = np.linalg.norm((atom + i) - (atom + f))
						if separation < (SiRadius + HRadius):
							dontAdd = True
					if not dontAdd:
						floatingAtoms.append([atom + f])
						added = True

					iterations += 1
					if iterations >= 20:
						added=True
						floatingBonds -= 1
				#f = f*SiHbondLength/np.linalg.norm(f)

			#Strained bonds
			elif numNeighbors == 4:
				strained = False
				for n, i in enumerate(bondVectors):
					strainedNeighbors = 0
					for j in bondVectors[:n] + bondVectors[n+1:]:
						angle = np.arccos(np.dot(i,j)/(np.linalg.norm(i)*np.linalg.norm(j)))
						bondangles.append(angle)
						angle_degrees = angle*180/np.pi
						if angle_degrees < 86 or angle_degrees > 131:
							#print("angle in degrees: ",angle_degrees)
							strainedNeighbors += 1
					if strainedNeighbors >= 2:
						#i is our highly strained bond, append H to it
						t = np.arccos(np.linalg.norm(i)/(2*SiHbondLength))
						H_alongbond = i*SiHbondLength/np.linalg.norm(i)
						rot_matrix = np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
						s = np.dot(rot_matrix, H_alongbond)
						
						strainedBonds += 1
						strained = True
						strainedAtoms.append([atom+s])

				if not strained:
					#non-defected == interstitial H
					if phase == "amorphous":
						interstitialSites += 1
						#determine whether it will be a dihydride or a single hydride
						r = np.random.random_sample()
						if r < percentDihydride:
							mono = False
						else:
							mono = True

						if mono:
							#simple method: negative of a Si bond
							i = np.random.randint(0,4)
							bond = -bondVectors[i]*SiHbondLength/np.linalg.norm(bondVectors[i])
							interstitialAtoms.append([atom + bond])
						else:
							#working solution: negative of two Si bonds
							i = np.random.randint(0,4)
							j = np.random.randint(0,4)
							while i == j:
								j = np.random.randint(0,4)
							bond1 = -bondVectors[i]*SiHbondLength/np.linalg.norm(bondVectors[i])
							bond2 = -bondVectors[j]*SiHbondLength/np.linalg.norm(bondVectors[j])
							interstitialAtoms.append([atom + bond1, atom+bond2])
					else:
						#only add H (for now) to the interstitial tetrahedral site
						compare = np.array([-0.5,-0.5,-0.5])
						compare = compare/np.linalg.norm(compare)
						#print("comparison: ",compare)
						for i in bondVectors:
							#print("vector: ",i/np.linalg.norm(i))
							if np.allclose(i/np.linalg.norm(i), compare,rtol = 1e-1):
								interstitialSites += 1
								inter = np.array([0.5*siUnitCellLength,0,0])
								interstitialAtoms.append([atom + inter])
								break

	else:
		for n, bondVectors in bondMap.items():
			atom = atomsXYZ[n]
			#Make the bond length of Si-H have a Gaussian width (to increase diversity of atomic neighborhoods)
			SiHbondLength = np.random.normal(SiHbondLength0,Sigma_bl)

			interstitialSites += 1
			numBonds = len(bondVectors)
			print("Number of bonds is: ",numBonds)
			if numBonds == 0:
				#determine whether it will be a dihydride or a single hydride
				r = np.random.random_sample()
				if r < percentDihydride:
					mono = False
				else:
					mono = True
				
				if mono:
					bondDirection = np.array([0.5,0.5,0.5])
					bond = SiHbondLength*bondDirection/np.linalg.norm(bondDirection)
					interstitialAtoms.append([atom + bond])
				else:
					bd1 = np.array([0.5,0.5,0])
					bd2 = np.array([0,0,0.5])
					b1 = SiHbondLength*bd1/np.linalg.norm(bd1)
					b2 = SiHbondLength*bd2/np.linalg.norm(bd2)
					interstitialAtoms.append([atom + b1, atom + b2])
			elif numBonds == 3:
				#Get the dangling bond vector
				a1 = bondVectors[0]
				a2 = bondVectors[1]
				a3 = bondVectors[2]
				d = np.cross(a1,a2) + np.cross(a2,a3) + np.cross(a3,a1)
				d = SiHbondLength*d/np.linalg.norm(d)
				interstitialAtoms.append([atom + d])
			elif numBonds >=5:
				#First try to add the H opposite to one of the floating bonds
				added = False
				for j in bondVectors:
					if not added:
						f = -SiHbondLength*j/np.linalg.norm(j)
						
						dontAdd = False
						for i in bondVectors:
							separation = np.linalg.norm((atom + i) - (atom + f))
							if separation < (SiRadius + HRadius):
								dontAdd = True
						if not dontAdd:
							interstitialAtoms.append([atom + f])
							added = True
				
				iterations = 0
				#Next, add it at random. If that fails, don't add it
				while not added:
					#Here there is no good spot to add H, except at random next to the central Si
					f0 = np.array([np.random.randint(-1,2),np.random.randint(-1,2),np.random.randint(-1,2)])
					fhat = f0/np.linalg.norm(f0)
					f = SiHbondLength*fhat
					dontAdd = False
					for i in bondVectors:
						separation = np.linalg.norm((atom + i) - (atom + f))
						if separation < (SiRadius + HRadius):
							dontAdd = True
					if not dontAdd:
						interstitialAtoms.append([atom + f])
						added = True

					iterations += 1
					if iterations >= 20:
						added=True
						interstitialSites -= 1
			else:
				#determine whether it will be a dihydride or a single hydride
				r = np.random.random_sample()
				if r < percentDihydride:
					mono = False
				else:
					mono = True

				if mono or numBonds == 1:
					#simple method: negative of a Si bond
					i = np.random.randint(0,numBonds)
					bond = -bondVectors[i]*SiHbondLength/np.linalg.norm(bondVectors[i])
					interstitialAtoms.append([atom + bond])
				else:
					#working solution: negative of two Si bonds
					#print("Added two interstitial H's!")
					i = np.random.randint(0,numBonds)
					j = np.random.randint(0,numBonds)
					while i == j:
						j = np.random.randint(0,numBonds)
					bond1 = -bondVectors[i]*SiHbondLength/np.linalg.norm(bondVectors[i])
					bond2 = -bondVectors[j]*SiHbondLength/np.linalg.norm(bondVectors[j])
					interstitialAtoms.append([atom + bond1, atom+bond2])

	#print("Bond angle distribution: ", np.mean(bondangles)*180/np.pi, " +- ", np.std(bondangles)*180/np.pi)
	
	#print("Strained bonds: ",strainedBonds)
	#print("Dangling bonds: ",danglingBonds)
	#print("Floating bonds: ",floatingBonds)
	#allNonDanglingDefects = floatingAtoms + strainedAtoms
	#allDefects = danglingAtoms + allNonDanglingDefects
	numDefects = strainedBonds + floatingBonds + danglingBonds


	#####################
	#ORDER OF IMPORTANCE
	#Dangling bonds
	#SiH2
	#Strained/ floating bonds
	#interstitial
	#####################

	#Gather all of the defected sites together
	nonInterstitialSites = dihydrideAtoms + danglingAtoms + strainedAtoms + floatingAtoms
	
	siteIndexes = []
	#Prioritize dihydride sites, then dangling bonds
	for k in range(dihydrideSites):
		siteIndexes.append(k)

	#Offset the danglingBonds indices by the number of dihydride sites
	danglingIndices = (np.random.permutation(danglingBonds) + dihydrideSites).tolist()
	siteIndexes.extend(danglingIndices)

	#Get the strained and floating bond sites
	strainedFloating = np.random.permutation(strainedBonds + floatingBonds)
	#Offset the strainedFloating indices by the number of dangling and dihydride sites
	strainedFloating = strainedFloating + danglingBonds + dihydrideSites
	siteIndexes.extend(strainedFloating.tolist())


	print("Total number of defects: ",numDefects, ". Composition is: 1) Dihydride sites- ",dihydrideSites,", 2) Dangling Bonds- ",danglingBonds, ", 3) Strained bonds- ",strainedBonds, ", 4) Floating bonds- ",floatingBonds)

	#Now decide how much H will be interstitial H
	#We will only add enough H to passivate either the total number of defects or to fill the desired atomic concentration, whichever is smaller
	#Liquid phase
	if phase == "liquid":
		num_H_add = num_H
		num_interstitial = num_H
	
	#Amorphous phase
	elif phase == "amorphous":
		num_H_add = min(num_H, numDefects)
		num_interstitial = math.floor((num_H_add - danglingBonds - 2*dihydrideSites)*percentInterstitial)
		if num_interstitial < 0:
			num_interstitial = 0
	
	#Diamond phases
	else:
		num_H_add = min(num_H, numDefects)
		num_interstitial = 0

		#In 50% of cases, include a random interstitial H no matter how much H has been added
		r = np.random.random_sample()
		if r <= 0.5:
			num_H_add += 1
			num_interstitial += 1

		#Add more interstitial H in the usual manner
		#addInterstitial = floor((num_H_add - danglingBonds - dihydrideSites)*percentInterstitial)
		num_interstitial += math.floor((num_H_add -num_interstitial- danglingBonds - 2*dihydrideSites)*percentInterstitial)
		if num_interstitial < 0:
			num_interstitial = 0

		print("Number of total H to add: ",num_H_add)
		print("Number of defects: ",numDefects)
		print("Number of dangling bonds: ", danglingBonds)
		print("Dihydride sites: ",dihydrideSites)
		print("Floating bonds: ",floatingBonds)
		print("Strained bonds: ",strainedBonds)
		print("Total amount of interstitial H to add: ",num_interstitial)

	#Get random permutation of interstitial indices
	interstitial_indexes = np.random.permutation(interstitialSites)

	H_atoms = []
	end = num_H_add - 1
	for i in range(num_H_add):
		if i <= end:
			if i < num_H_add - num_interstitial:
				item = nonInterstitialSites[siteIndexes[i]]
			else:
				item = interstitialAtoms[interstitial_indexes[i - (num_H_add - num_interstitial)]]

			if len(item) == 1:
				H_atoms.append(item[0])
			else:
				if i != end:
					H_atoms.append(item[0])
					H_atoms.append(item[1])

					#Added two H's so one less will be added total
					end -= 1

				#Only one more H can be added, so prevent a dihydride addition
				else:
					j = i
					while len(item) != 1:
						j += 1
						if i < num_H_add - num_interstitial:
							item = nonInterstitialSites[siteIndexes[j]]
						else:
							item = interstitialAtoms[interstitial_indexes[j - (num_H_add - num_interstitial)]]
					H_atoms.append(item[0])

	return H_atoms

#################################
# Now we actually run the code! #
#################################

#Get the input filename
inputFileName = sys.argv[1]
outputFileTree = inputFileName[:-4]

#First get the file containing all the structures
directory = os.getcwd()
file_extension = '/Pure Si structures/'
inputfile = directory + file_extension + inputFileName

#Now start to read in the LAMMPS file, and create the &SYSTEM namelist
lines = open(inputfile,"r").read().splitlines()
line_iterator = 0
file_iterator = 1
nAtoms = 1
infoDict = {}
atomsExtended = []
atomsXYZ = []
for n,line in enumerate(lines):
	#get the number of atoms in the structure
	if line_iterator == 0:
		nAtoms = int(line)
	#Read in the header information, and put it into a dictionary
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
					if n%3 == 0:
						item = l
					elif n%3 == 2:
						end = start + int(l)
						fielddata[item] = [start,end]
						start = end
				infoDict[fieldname] = fielddata
			else: 
				infoDict[fieldname] = fields[1]
	#Read in the extended xyz atom data
	elif line_iterator <= nAtoms + 1:
		atomInfo = line.split()
		atomsExtended.append(atomInfo)
		xyz_start = infoDict['Properties']['pos'][0]
		xyz_end = infoDict['Properties']['pos'][1]
		atomsXYZ.append([float(s) for s in atomInfo[xyz_start:xyz_end]])
	line_iterator += 1

	if line_iterator > nAtoms + 1:
		print("File number: ",file_iterator)

		#We have reached the end of the structure, so now we need to do 2 things:
		#First, we need to create our map of bond lengths
		#Second, we need to add a specified number of H to our Si structure
		#Third, we need to create the DFT input file
		#Let's do the first task:
		latticeComponents = infoDict['Lattice'].split()
		latticeVectors = np.reshape(latticeComponents,(3,3))
		bondMap = getDistanceMap(atomsXYZ, latticeVectors, SiSibondLength)

		#Get k-point numbers
		a1 = np.asarray([float(i) for i in latticeVectors[0]])
		a2 = np.asarray([float(i) for i in latticeVectors[1]])
		a3 = np.asarray([float(i) for i in latticeVectors[2]])

		a1Norm = np.linalg.norm(a1)
		a2Norm = np.linalg.norm(a2)
		a3Norm = np.linalg.norm(a3)

		k1 = math.ceil(2*math.pi/(a1Norm*0.2))
		k2 = math.ceil(2*math.pi/(a2Norm*0.2))
		k3 = math.ceil(2*math.pi/(a3Norm*0.2))

		#Now add H to our Si structure
		num_H = (int)(nAtoms*percent_H/(1-percent_H))
		atomsXYZ_H = getH_XYZ(atomsXYZ, bondMap, num_H)
		natoms_total = nAtoms + len(atomsXYZ_H)

		#Now we need to create the DFT input file
		#Create the output file
		outputDirectory = directory + "/" + "Si:H_structure_QE_input_files/" + outputFileTree + "/"
		if not os.path.exists(outputDirectory):
			os.makedirs(outputDirectory)
		outputfile = outputDirectory + phase + str(file_iterator) + ".in"

		#Now create the &CONTROL namelist
		out = open(outputfile,'w')
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
		out.write("  ibrav = 0,")
		#out.write(" celldm(1) = " + "{:.9f}".format(xSize*1.889716164) + ",")
		#out.write(" celldm(3) = " + "{:.9f}".format(zRatio) + ",")
		out.write(" nat = " + str(natoms_total) + ",")
		out.write(" ntyp = 2,\n")
		out.write("  ecutwfc\t= " + ecutwfc + ",\n")
		#out.write("  nbnd\t= 1728,\n")
		out.write("  degauss\t= " + gauss +",\n")
		out.write("  occupations\t= " + occupation + ",\n")
		out.write("  smearing\t= " + smearing + ",\n")
		out.write("/\n")

		#Now create the &ELECTRONS namelist
		out.write("&ELECTRONS\n")
		out.write("  electron_maxstep = 200,\n")
		out.write("  conv_thr = 1.D-8,\n")
		out.write("  mixing_beta = 0.3D0,\n")
		out.write("  scf_must_converge=.false.,\n")
		out.write("/\n")

		#Now create the &IONS namelist
		out.write("&IONS\n")
		out.write("/\n")

		#Now create the ATOMIC_SPECIES card
		out.write("ATOMIC_SPECIES\n")
		out.write("Si  28.08  " + pseudo_potential_Si + "\n")
		out.write("H  1.00784  " + pseudo_potential_H + "\n")

		#Now create the CELL_PARAMETERS card
		out.write("CELL_PARAMETERS { angstrom }\n")
		for i in latticeVectors:
			out.write("   " + i[0] + "  " + i[1] + "  " + i[2] + "\n")

		out.write("K_POINTS {automatic}\n")
		out.write(str(k1) + " " + str(k2) + " " + str(k3) + "  0 0 0\n")

		out.write("ATOMIC_POSITIONS (angstrom)\n")
		for i in atomsXYZ:
			out.write("Si "+ "{:.10f}".format(i[0]) + " " + "{:.10f}".format(i[1]) + " " + "{:.10f}".format(i[2]) + "\n")

		for i in atomsXYZ_H:
			out.write("H "+ "{:.10f}".format(i[0]) + " " + "{:.10f}".format(i[1]) + " " + "{:.10f}".format(i[2]) + "\n")

		#reset the iterator
		line_iterator = 0
		file_iterator += 1
		atomsXYZ = []
