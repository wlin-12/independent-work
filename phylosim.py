import os
from ete3 import Tree
from TreeUtils import findDomains
import sys
import Similarity
import glob
import numpy as np
import string
import random
import matplotlib.pyplot as plt

def findAndAlign(hmmfile, sequence):
	possStarts, possEnds, possSequences = findDomains(sequence, hmmfile)
	starts = []
	ends = []
	sequences = []
	for i in range(len(possSequences)):
		if len(possSequences[i]) == 23 and any(base.islower() for base in possSequences[i]) == False:
			starts.append(possStarts[i])
			ends.append(possEnds[i])
			sequences.append(possSequences[i])

	origCount = len(sequences)
	# print("ORIG: " + str(origCount))

def uniqueID(fileName):
	fasta = fastaToDictionary(fileName)
	ofile = open("data/U_" + fileName, "w")
	for entry in fasta:
		ofile.write(">" + entry + "_" + fileName + "\n" + fasta[entry] + "\n")
	ofile.close()

def remove(fileName, unrooted, rooted, data, branchLengths):
	unrootedFiles = ['result.', 'parsimony.', 'log.', 'info.', 'bestTree.', 'bestTree.R']
	rootedFiles = ['nodeLabelledRootedTree.', 'marginalAncestralStates.', 'marginalAncestralProbabilities.', 'info.']

	if unrooted == True:
		for i in unrootedFiles:
			files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_unrooted/RAxML_' + i + fileName)
			for f in files:
				os.remove(f)
	if rooted == True:
		for i in rootedFiles:
			files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_rooted/RAxML_' + i + fileName)
			for f in files:
				os.remove(f)
	if data == True:
		files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/data/OUTGROUP_' + fileName)
		for f in files:
			os.remove(f)
		files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/data/U_' + fileName)
		for f in files:
			os.remove(f)
	if branchLengths == True:
		files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_BL' + fileName)
		for f in files:
			os.remove(f)


def fastaToDictionary(fileName):
	fasta = {}
	with open('data/' + fileName) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				active_sequence_name = line[1:]
				if active_sequence_name not in fasta:
					fasta[active_sequence_name] = []
				continue
			sequence = line
			fasta[active_sequence_name].append(sequence)
	for entry in fasta:
		value = ""
		for line in fasta[entry]:
			value = value + line
		fasta[entry] = value
	return fasta

def generateRootedTree(sequenceFileName):
	remove(sequenceFileName, True, True, True, False)
	fasta = fastaToDictionary(sequenceFileName)

	entries = []
	for entry in fasta:
		entries.append(entry)
	ofile = open("data/OUTGROUP_" + sequenceFileName, "w")
	length = 0
	for entry in fasta:
		ofile.write(">" + entry + "\n" + fasta[entry] + "\n")
		length = len(fasta[entry])
	OUTGROUP = "A" * length
	ofile.write(">" + 'OUTGROUP' + "\n" + OUTGROUP + "\n")
	ofile.close()

	# Build unrooted tree containing fictional outgorup
	relativeURL = 'data/' + sequenceFileName
	outputFileName = sequenceFileName
	substitutionModel = "PROTGAMMABLOSUM62"
	outputDirectory = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_unrooted/"
	threads = "10"
	randomSeed = np.random.randint(0,2000)
	command = "./raxml -s data/" + sequenceFileName + " -n " + outputFileName + " -m " + substitutionModel + " -w " + outputDirectory + " -T " + threads + " -p " + str(randomSeed) + " 1>>log_file"
	os.system(command)

	#root based on midpoint
	current = Tree('trees_unrooted/RAxML_bestTree.'  + sequenceFileName)
	current.set_outgroup(current.get_midpoint_outgroup())
	current.prune(entries)

	ofile = open('trees_unrooted/RAxML_bestTree.R_'  + sequenceFileName , "w")
	ofile.write(current.write(format = 1))
	ofile.close()

def generateMarginalAncestralStates(sequenceFileName):
	remove(sequenceFileName, True, True, False, False)
	rootedTree = 'trees_unrooted/RAxML_bestTree.R_'  + sequenceFileName
	randomSeed = np.random.randint(0,2000)
	outputFileName = sequenceFileName
	substitutionModel = "PROTGAMMABLOSUM62"
	threads ="10"
	outputDirectory = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_rooted/"
	command = "./raxml -f A -t " + rootedTree + " -s data/" + sequenceFileName + " -n " + outputFileName + " -m " + substitutionModel + " -w " + outputDirectory + " -T " + threads + " -p " + str(randomSeed) + " 1>>log_file"
	os.system(command)

	return outputDirectory + "RAxML_marginalAncestralStates." + sequenceFileName

def readMarginalAncestralStates(inputFile):
	internalNodes = {}
	with open(inputFile, 'r') as file_one:
		for line in file_one:
			data = line.split(' ')
			internalNodes[data[0]] = data[1].strip()
	return internalNodes

def sequenceSimilarities(subtrees):
	# print(subtrees)
	length = len(subtrees)
	maxSimilarity = -1
	maxTreeX = ""
	maxNodeX = ""
	maxTreeY = ""
	maxNodeY = ""
	for fileA, sequencesA in sorted(subtrees.iteritems()):
		for fileB, sequencesB in sorted(subtrees.iteritems()):
			# print(fileA)
			# print(fileB)
			if fileA != fileB:
				root = subtrees[fileA]['ROOT']
				for sequence in subtrees[fileB]:
					similarity = Similarity.domainSim(root,subtrees[fileB][sequence])
					if similarity > maxSimilarity:
						maxSimilarity = similarity
						maxTreeX = fileA
						maxNodeX = 'ROOT'
						maxTreeY = fileB
						maxNodeY = sequence
	# indices = sorted([(maxTreeX, maxNodeX), (maxTreeY, maxNodeY)], key=lambda x: x[0])
	return maxTreeX, maxTreeY, maxNodeY, maxSimilarity

def randomString(stringLength=8):
	letters = string.ascii_lowercase
	return ''.join(random.choice(letters) for i in range(stringLength))

def joinFasta(file1, file2):
	joinedFastaName = file1.split('.')[0] + "-" + file2
	data = data2 = "" 
	with open('data/' + file1) as fp: 
		data = fp.read() 
	with open('data/' + file2) as fp: 
		data2 = fp.read() 
	data += "\n"
	data += data2 
	with open ('data/' + joinedFastaName, 'w') as fp: 
		fp.write(data)

def joinSubtree(maxTreeX, maxTreeY, maxNodeY, maxSimilarity, fastaFiles, connectMethod):
	file1Tree = Tree('trees_rooted/RAxML_nodeLabelledRootedTree.'  + maxTreeX, format = 8)
	file2Tree = Tree('trees_rooted/RAxML_nodeLabelledRootedTree.'  + maxTreeY, format = 8)

	# print("MOST RECENT JOIN: " + maxTreeX + " and " + maxTreeY)
	# print 
	joinedFastaName = maxTreeX.split('.')[0] + "-" + maxTreeY
	# print(file2Tree.write(format=8))
	# print(maxNodeY)
	if connectMethod == 0:
		target = file2Tree.search_nodes(name=maxNodeY)
		if target[0].is_root():
			name =randomString(5)
			root = Tree(name + ";")
			root.add_child(file1Tree)
			root.add_child(file2Tree)
			ofile = open('trees_unrooted/RAxML_bestTree.R_'  +  joinedFastaName, "w")
			ofile.write(root.write(format = 1))
			ofile.close()
		else:
			parent = target[0].up
			target[0].detach()
			new = parent.add_child(name=randomString(5))
			new.add_child(file1Tree)
			new.add_child(target[0])		
			ofile = open('trees_unrooted/RAxML_bestTree.R_'  +  joinedFastaName, "w")
			ofile.write(file2Tree.write(format = 1))
			ofile.close()
	# if connectMethod == 1:
	# 	name =randomString(5)
	# 	root = Tree(name + ";")
	# 	root.add_child(file1Tree)
	# 	root.add_child(file2Tree)

	# 	ofile = open('trees_unrooted/RAxML_bestTree.R_'  +  joinedFastaName, "w")
	# 	print(root.write(format=1))
	# 	ofile.write(root.write(format = 1))
	# 	ofile.close()

	fastaFiles.remove(maxTreeX)
	fastaFiles.remove(maxTreeY)
	fastaFiles.append(joinedFastaName)
	joinFasta(maxTreeX, maxTreeY)

	return fastaFiles

def joinTree(fastaFiles, connectMethod):
	allStates = {}
	for file in fastaFiles:
		# print(file)
		output = generateMarginalAncestralStates(file)
		# print("LADIDA")
		allStates[file] = readMarginalAncestralStates(output)
	# print(allStates)
	maxTreeX, maxTreeY, maxNodeY, maxSimilarity = sequenceSimilarities(allStates)
	return joinSubtree(maxTreeX, maxTreeY, maxNodeY, maxSimilarity, fastaFiles, connectMethod)

def generateTree(fastaFiles, connectMethod):
	while len(fastaFiles) > 1:	
		# print("Current iteration:")
		# print(fastaFiles)
		remainingFiles = joinTree(fastaFiles, connectMethod)
		return generateTree(remainingFiles, connectMethod)
	generateMarginalAncestralStates(fastaFiles[0])
	return fastaFiles[0]

def branchLengths(sequenceFileName):
	rootedTree = 'trees_unrooted/RAxML_bestTree.R_'  + sequenceFileName
	outputFileName = sequenceFileName
	substitutionModel = "PROTGAMMABLOSUM62"
	threads ="10"
	outputDirectory = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_BL/"
	command = "./raxml -f e -t " + rootedTree + " -s data/" + sequenceFileName + " -n " + outputFileName + " -m " + substitutionModel + " -w " + outputDirectory + " -T " + threads + " 1>>log_file"
	os.system(command)

	return outputDirectory + "RAxML_marginalAncestralStates." + sequenceFileName

def main():
	# files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_BL/*')
	# for f in files:
	# 	os.remove(f)
	# files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_unrooted/*')
	# for f in files:
	# 	os.remove(f)
	# files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_rooted/*')
	# for f in files:
	# 		os.remove(f)
	# for i in range(0,1):
	# 	files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/data/OUTGROUP_*')
	# 	for f in files:
	# 		os.remove(f)
	# 	files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/log_file')
	# 	for f in files:
	# 		os.remove(f)
	# 	directory = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/data"
	# 	data = []
	# 	for filename in os.listdir(directory):
	# 		# print(filename)
	# 		# if filename.endswith("_" + str(i)): 
	# 		data.append(filename)
	# 	# print(data)
	# 	# print(data)
	# 	#data = ['0.fasta', '1.fasta', '2.fasta', '3.fasta','4.fasta']
	# 	#identified = []
	# 	# for file in data:
	# 	# 	uniqueID(file)
	# 	# 	identified.append("U_" + file)
	# 	for file in data:
	# 		generateRootedTree(file)
		# print("AAAAAAAAA")
		# nodeLabelledRootedTree = generateTree(data,0)
		# # remove(nodeLabelledRootedTree, False, False, False, True)
		# branchLengthTree = branchLengths(nodeLabelledRootedTree)


	output = {}
	directory1 = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_BL/"
	for i in range(0,150):
		# if i not in [28,108,66,285,388,432,472,492]:
			# print(filename)
		for fileName in os.listdir(directory1):
			if fileName.endswith("_" + str(i) + ".fasta"): 
				tag = fileName.split(".")[0]
				joined = fileName.split(".")[1] + "." + fileName.split(".")[2]
				output[i] = joined
			# parts = filename.split("_")
			# if len(parts) == 4 and filename.endswith("_" + str(i) + ".fasta"): 
			# 	output[i] = parts[2] + "_" + parts[3]
	# print(output)
	directory2 = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/trees_rooted/RAxML_nodeLabelledRootedTree."
	directory3 = "/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_guestTrees/"

	# print(output)
	rf_values = []
	compare = []
	for step in output:
		# print(output[step])
		print(directory2 + output[step])
		t1 = Tree(directory2 + output[step], format = 1)
		t2 = Tree(directory3 + "guestTrees_step_" + str(step), format = 1)
		rf = t1.robinson_foulds(t2)
		rfN = t1.compare(t2)["norm_rf"]
	  	rf_values.append(rf[0])
	  	compare.append(rfN)
	
	# print(rf_values)
	# print(compare)

	# print("----------------------------------------------------------------------")

	print("Average Distance: " + str(np.average(rf_values)))
	print("Variance of Distance " + str(np.var(rf_values)))
	print("Average Distance: " + str(np.average(compare)))
	print("Variance of Distance " + str(np.var(compare)))


	plt.hist(rf_values, bins=[0, 20, 40, 60, 80, 100, 120,140])
	plt.title('Distribution of Robinson-Foulds Distances')
	plt.xlabel('Robinson-Foulds Distances')
	plt.ylabel('Number of Reconstructed Trees')
	plt.xlim((0, 150))
	plt.ylim((0, 60))
	plt.show()

	plt.hist(compare, bins=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
	plt.title('Distribution of Normalized Robinson-Foulds Distances')
	plt.xlabel('Robinson-Foulds Distances')
	plt.ylabel('Number of Reconstructed Trees')
	plt.show()


if __name__ == "__main__":
	main()


