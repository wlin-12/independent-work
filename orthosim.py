from ete3 import Tree
from rootSequence import genRandomSequence2
from TreeUtils import findDomains
from Utils import sortBy
from CreateSequences import evolveLinker
import numpy as np
from math import log
import random
import HostTreeGen, GuestTreeGen
from Utils import printProgressBar
from TreeUtils import s3
from stats import expfunc, gaussNoise
import numpy as np
import Similarity
import matplotlib.pyplot as plt
import RetrieveOrthogroup
import glob
import os
import string

host = Tree("(C,(A,B)D)E;", format = 8)
guest = Tree("(((a:0.2,r:0.2)p:0.3,b:0.5)e:0.5,((c:0.4,s:0.4)q:0.3,d:0.1)f:0.3)z;", format = 1)

nodemap = {host&"A": [guest&'a',guest&'s'],host&'B': [guest&'c',guest&'r'],host&'C': [guest&'b',guest&'d'],host&'D':[guest&'p',guest&'q'],host&'E': [guest&'e',guest&'f',guest&'z']}
print("HOST TREE (W/ NAMES):")
print(host.get_ascii(attributes=['name']))
print
print("GUEST TREE (W/ NAMES):")
print(guest.get_ascii(attributes=['name']))
print
print("GUEST TREE (W/ DISTANCES):")
print(guest.get_ascii(attributes=['dist']))
print

for node in guest.traverse():
    if node.name in ['z']:
        node.add_feature('event', "DUPLICATION")
    elif node.name in ['d']:
        node.add_feature('event', "LOSS")
    else:
        node.add_feature('event', "SPECIATION")

print("GUEST TREE (W/ EVENTS):")
print(guest.get_ascii(attributes=['event']))
print
domains = 11
sequence = genRandomSequence2(domains)

def findLinkers(starts, ends, startingSequence):
    linkerStarts = []
    linkerEnds = []
    linkerSequences = []
    for i in range(len(starts)):
        if i == 0:
            if starts[0] != 0:
                linkerStarts.append(0)
                linkerEnds.append(starts[0]-1)
                linkerSequences.append(startingSequence[0:starts[0]])
        else:
            linkerStarts.append(ends[i-1]+1)
            linkerEnds.append(starts[i]-1)
            linkerSequences.append(startingSequence[ends[i-1]+1:starts[i]])
    if ends[len(ends)-1] != len(startingSequence)-1:
        start_position = ends[len(ends)-1] + 1
        end_position = len(startingSequence)
        linkerStarts.append(start_position)
        linkerEnds.append(end_position)
        linkerSequences.append(startingSequence[start_position:])
    return linkerStarts, linkerEnds, linkerSequences

def reconstructSequence(starts, ends, sequences, linkerStarts, linkerEnds, linkerSequences):
    sequence = ""
    if starts[0] < linkerStarts[0]:
        if len(starts) > len(linkerStarts):
            for i in range(len(starts)-1):
                sequence += sequences[i] + linkerSequences[i]
            sequence += sequences[len(starts)-1]
        else:
            for i in range(len(starts)):
                sequence += sequences[i] + linkerSequences[i]
    else:
        if len(linkerStarts) > len(starts):
            for i in range(len(linkerStarts)-1):
                sequence += linkerSequences[i] + sequences[i]
            sequence += linkerSequences[len(linkerSequences)-1]
        else:
            for i in range(len(linkerStarts)):
                sequence += linkerSequences[i] + sequences[i]
    return sequence

def probability(score):
     return 0.01516616359*2**(score/1000)

def mutateDomain(domain,hmmfile,branchlength):
    positions = 23
    f = open(hmmfile, "r")
    hmmPositions = []
    for line in f:
        hmmPositions.append(line)
    hmmPositions= hmmPositions[20:len(hmmPositions)-2]
    parsed = []
    for i in range(len(hmmPositions)):
        if (i+1)%3 == 1:
            position = hmmPositions[i]
            scores = position.split()
            parsed.append(scores[1:-1])

    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    probabilityDistribution = []
    for i in range(positions):
        for j in range(len(aminoacids)):
            parsed[i][j] = probability(int(parsed[i][j]))
        scale = 1 / np.sum(parsed[i])
        for j in range(len(aminoacids)):
            parsed[i][j] *= scale
        probabilityDistribution.append(parsed[i])

    cumulativeProbabilityDistribution = []
    for i in range(positions):
        cumSum = 0
        current = []
        for j in range(len(aminoacids)):
            cumSum += probabilityDistribution[i][j]
            current.append(cumSum)
        cumulativeProbabilityDistribution.append(current)
    # print(cumulativeProbabilityDistribution)

    informationEntropy = []
    for i in range(len(probabilityDistribution)):
        sum = 0
        for j in range(len(probabilityDistribution[i])):
            sum += -probabilityDistribution[i][j] * log(probabilityDistribution[i][j],2)
        informationEntropy.append((i,sum))
    informationEntropy = sorted(informationEntropy, key=lambda x: x[1])

    # print(domain)
    mutations = round(branchlength * positions)
    for i in range(int(mutations)):
        position = informationEntropy[i][0]
        current = cumulativeProbabilityDistribution[position]
        value = random.uniform(0,.999)
        for j in range(len(cumulativeProbabilityDistribution[i])):
            if value < cumulativeProbabilityDistribution[i][j]:
                domain = domain[0:j] + aminoacids[j] + domain[j+1:]
                # print(domain)
                break
    return domain

def removeSimulations():
    files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/data/*')
    for f in files:
        os.remove(f)
    files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_guestTrees/*')
    for f in files:
        os.remove(f)
    files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_marginalAncestralStates/*')
    for f in files:
        os.remove(f)
    files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_distances/*')
    for f in files:
        os.remove(f)
    files = glob.glob('/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_leaves/*')
    for f in files:
        os.remove(f)

def domainEvolution(host, guest, hmmfile, nodemap, sequence):
    possStarts, possEnds, possSequences = findDomains(sequence, hmmfile)
    starts = []
    ends = []
    sequences = []
    positions = {}
    for i in range(len(possSequences)):
        if len(possSequences[i]) == 23 and any(base.islower() for base in possSequences[i]) == False:
        # if len(possSequences[i]) == 23 and any(base.islower() for base in possSequences[i]) == False:
            starts.append(possStarts[i])
            ends.append(possEnds[i])
            sequences.append(possSequences[i])

    origCount = len(sequences)
    print("ORIG: " + str(origCount))
    linkerStarts, linkerEnds, linkerSequences = findLinkers(starts, ends, sequence)

    orthogroup = []
    bookkeeping = {}
    internalNodes = {}

    guestData = {}
    print("*******************************************")
    print(host.get_tree_root())
    for node in host.traverse():
        print("HOST NODE: " + node.name)
        print("HOST NODE CHILDREN: " + str(node.children))

        # find the root sequence
        mapped_guest = nodemap[node]
        print("NODEMAP: " + str(mapped_guest))
        guestRoots = []

        if node.is_root():
            guestRoots.append(guest.get_tree_root())
        else:
            for guestNode in mapped_guest:

                if guestNode.up not in mapped_guest:
                    guestRoots.append(guestNode)

        print("GUEST ROOTS: " + str(guestRoots))

        # # map guestTree node to domains
        if node.is_root():
            node.add_feature('sequences', sequences[:])
            node.add_feature('starts', starts[:])
            node.add_feature('ends', ends[:])
            node.add_feature('positions', dict(positions))
            node.add_feature('linkerStarts', linkerStarts[:])
            node.add_feature('linkerEnds', linkerEnds[:])
            node.add_feature('linkerSequences', linkerSequences[:])
            for i in range (len(guestRoots)):
                guestRoots[i].add_feature('sequences', sequences[:][i])
                guestRoots[i].add_feature('starts', starts[:][i])
                guestRoots[i].add_feature('ends', ends[:][i])
        else:
            betweenNodeLinkers = node.up.linkerSequences[:]
            betweenNodeSequences = node.up.sequences[:]
            distance = node.up.get_distance(node)
            # sequence level linker evolution
            for i in range(len(betweenNodeLinkers)):
                betweenNodeLinkers[i] = evolveLinker(betweenNodeLinkers[i], distance)

            # domain level evolution
            for i in range(len(betweenNodeSequences)):
                betweenNodeSequences[i] = mutateDomain(betweenNodeSequences[i], hmmfile, distance)
            node.add_feature('sequences', betweenNodeSequences[:])
            node.add_feature('starts', node.up.starts[:])
            node.add_feature('ends', node.up.ends[:])
            node.add_feature('positions',dict(node.up.positions))
            node.add_feature('linkerStarts', node.up.linkerStarts[:])
            node.add_feature('linkerEnds', node.up.linkerEnds[:])
            node.add_feature('linkerSequences', betweenNodeLinkers[:])

        current_sequences = node.sequences
        current_starts = node.starts
        current_ends = node.ends
        current_positions = node.positions
        current_linkerStarts = node.linkerStarts
        current_linkerEnds = node.linkerEnds
        current_linkerSequences = node.linkerSequences

        print("INPUT POSITIONS: " + str(current_positions))

        print("INPUT: " + str(sequences))

        print("MAPPING: " + str(mapped_guest))

        # initialize positions 
        pos_init = 0

        extra = dict({})

        print("-------------------------------------------")
        for root in guestRoots:
            print("CURRENTLY EXAMINING GUEST ROOT: " + root.name)
            # for a subtree by traversing from a root node,
            # see if node has chilren that are in the list,
            # cut off when it doesn't to form subtree

            subtree = root.copy("deepcopy")
            # print(subtree.write(format=8))
            for newbie in subtree.iter_descendants():
                if newbie not in mapped_guest:
                    newbie.detach()
            # print(subtree.write(format=8))
            # initialize distances
            distances = []
            distances.append([root,0])  
            closestNode = distances[0][0]
            closestDistance = distances[0][1]
            index = 0

            # obtain domain information (to be updated later)
            if node.is_root():
                root_sequences = root.sequences
                root_starts = root.starts
                root_ends = root.ends
                # index = guestRoots.index(root)
            else:
                for i in current_positions:
                    if i.name == root.up.name:
                        index = current_positions[i]
                        del current_positions[i]
                        current_positions[root] = index
                # index = current_positions[root.up]
                # if root.up in current_positions:
                #     del current_positions[root.up]
                root_sequences = current_sequences[index]
                root_starts = current_starts[index]
                root_ends = current_ends[index]

            if pos_init == 0 and node.is_root():
                current_positions[root] = 0
                pos_init = 1

            length = len(root_sequences)

            count = 1

            # iterate by minimum distance
            while True:
                print("EVENT TITLE: " + closestNode.event)
                # print("CURRENT POSITIONS: " + str(current_positions))
                # print(index)
                # check event node, update positions list
                if closestNode.event == "DUPLICATION":
                    for position in current_positions:
                        if current_positions[position] > index:
                            current_positions[position] += 1
                    if node.is_root and root.is_root:
                        oldPosition = current_positions[closestNode]
                        del current_positions[closestNode]
                    else:
                        oldPosition = current_positions[closestNode.up]
                        del current_positions[closestNode.up]
                    current_positions[closestNode.children[0]] = oldPosition
                    current_positions[closestNode.children[1]] = oldPosition + 1
                    linkerLength = 0
                    if index == 0:
                        linkerLength = 5
                    else:
                        linkedLength = len(current_linkerSequences[index])

                    current_starts.append(root_starts + linkerLength + length)
                    current_ends.append(root_ends + linkerLength + length)
                    current_sequences.append(root_sequences)
                    current_linkerStarts.append(current_linkerStarts[index] + linkerLength + length)
                    current_linkerEnds.append(current_linkerEnds[index] + linkerLength + length)
                    current_linkerSequences.append(current_linkerSequences[index])

                elif closestNode.event == "LOSS":
                    current_starts.pop(index)
                    current_ends.pop(index)
                    current_sequences.pop(index)
                    current_linkerStarts.pop(index)
                    current_linkerEnds.pop(index)
                    current_linkerSequences.pop(index)
                    for position in current_positions:
                        if current_positions[position] > index:
                            current_positions[position] -= 1
                    del current_positions[current_positions.keys()[index]]
                elif closestNode.event == "SPECIATION":
                    closestNodeSearch = guest&closestNode.name
                    if closestNodeSearch.up in current_positions:
                        del current_positions[closestNodeSearch.up]
                        current_positions[closestNode] = index

                # # sort to maintain original order in case of duplication
                current_sequences = sortBy(current_sequences, current_starts)
                current_starts.sort()
                current_ends.sort()
                current_positions = dict(sorted(current_positions.items(), key=lambda x: x[1]))
                current_linkerSequences = sortBy(current_linkerSequences, current_linkerStarts)
                current_linkerStarts.sort()
                current_linkerEnds.sort()

                print("POSITIONS: " + str(current_positions))
                print("STARTS: " + str(current_starts))
                print("ENDS: " + str(current_ends))

                print("CLOSEST NAME AND DISTANCE:")
                print(closestNode.name)
                print(closestDistance)
                guestNode = closestNode.name
                guestData[closestNode.name] = [current_sequences[index], current_starts[index], current_ends[index]]

                print("PRE-DELETION DISTANCES:")
                print(distances)

                del distances[distances.index([closestNode,closestDistance])]
                print("POST-DELETION DISTANCES:")
                print(distances)
                closestChildren = closestNode.children
                print("CLOSEST CHILDREN:")
                print(closestChildren)

                # subtract distance to closest from every remaining considered gene
                for gene in distances:
                    gene[1] -= closestDistance
                # add the closest's childrento the list
                for child in closestChildren:
                    # print(mapped_guest)
                    if child in mapped_guest:
                        distances.append([child,child.dist])
                print("UNSORTED DISTANCES: " + str(distances))
                sortedDistances = sorted(distances, key=lambda x: x[1])
                print("SORTED DISTANCES: " + str(sortedDistances))
                if len(sortedDistances) > 0:
                    closestNode = sortedDistances[0][0]
                    closestDistance = sortedDistances[0][1]
                    # sequence level linker evolution
                    for i in range(len(current_linkerSequences)):
                        current_linkerSequences[i] = evolveLinker(current_linkerSequences[i], closestDistance)

                    # domain level evolution
                    for i in range(len(current_sequences)):
                        old = current_sequences[i]
                        current_sequences[i] = mutateDomain(current_sequences[i], hmmfile, closestDistance)

                    distances = sortedDistances
                    sequence = reconstructSequence(current_starts, current_ends, current_sequences, current_linkerStarts, current_linkerEnds, current_linkerSequences)
                    internalNodes[guestNode] = sequence
                    guestData[closestNode.name] = [current_sequences[index], current_starts[index], current_ends[index]]

                    print("-------------------------------------------")
                    current_positions[closestNode] = index + 1
                    index += 1
                    if closestNode not in mapped_guest:
                        break
                else:
                    distances = sortedDistances
                    sequence = reconstructSequence(current_starts, current_ends, current_sequences, current_linkerStarts, current_linkerEnds, current_linkerSequences)
                    internalNodes[guestNode] = sequence
                    guestData[closestNode.name] = [current_sequences[index], current_starts[index], current_ends[index]]
                    print("-------------------------------------------")
                    break


            print("Finished examining guest root: " + root.name)

            # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        node.add_feature("positions", current_positions)
        node.add_feature("sequences", current_sequences)
        node.add_feature("starts", current_starts)
        node.add_feature("ends", current_ends)

        if node.is_leaf() == 1:
            bookkeeping[node.name] = [current_sequences, current_starts, current_ends]
        print("*******************************************")
        print("ENDING sequences: " + str(current_sequences))
        finalSequence = reconstructSequence(current_starts, current_ends, current_sequences, current_linkerStarts, current_linkerEnds, current_linkerSequences)
        node.add_feature("final", finalSequence)
        print
        print

        if node.is_leaf() == True:
            if len(finalSequence) > 10:
                orthogroup.append(finalSequence)
            # orthogroup[node.name] = finalSequence
        # print(orthogroup)


    return orthogroup, bookkeeping, internalNodes, guestData

def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def cumulative(tree):
    dist = {}
    for node in tree.traverse():
        if node.is_root() == False and node.name != "G0_0":
            distance = node.dist + dist[node.up.name]
            dist[node.name] = distance
        else:
            dist[node.name] = node.dist
    return dist
        

# output = domainEvolution(host, guest, 'zf-C2H2-232.hmm', nodemap, sequence)
# print("ORTHOGROUP IS: " + str(output))
# averages = []
# variances = []
# notDegenerate = 0
# removeSimulations()
for step in range(200,220):
    print("STEP: " + str(step))
    sizes = []
    treeHeight = 0.35
    while True:
        host = HostTreeGen.birthDeathTree(.3, .1, treeHeight)
        if len(host) > 25:
            break
    # print(host.write())
    guest, nodemap = GuestTreeGen.buildGuestTree(host, s3, expfunc, treeHeight / 5, gaussNoise, 10)

    # print(nodemap)
    test = []
    for node in guest.children[0].traverse("postorder"):
      # Do some analysis on node
      test.append(node)

    # print(test)
    # # print(nodemap)
    # for node in host.traverse():
    #     print(str(node) + " " + str(node.is_root()))

    prunedMap = {}
    for i in nodemap:
        if i not in prunedMap:
            prunedMap[i] = []
        for mapped in nodemap[i]:
            if mapped in test:
                prunedMap[i].append(mapped)
    print(prunedMap)

    # print(host.write(format=8))
    # print(nodemap)
    # print(guest.write(format=8))
    # print(guest.children)
    lucky = guest.children[0].copy("deepcopy")

    print("HOST TREE (W/ NAMES):")
    print(host.get_ascii(attributes=['name']))
    print
    print(host.get_tree_root())
    print
    print("GUEST TREE (W/ NAMES):")
    print(guest.children[0].get_ascii(attributes=['name']))
    print
    print("GUEST TREE (W/ DISTANCES):")
    print(guest.children[0].get_ascii(attributes=['dist']))
    print
    print("GUEST TREE (W/ EVENTS):")
    print(guest.children[0].get_ascii(attributes=['event']))
    print

    # print(lucky.write(format=8))
    outputOrthogroup, bookkeeping, internalNodes, guestData = domainEvolution(host, lucky, 'zf-C2H2-232.hmm', prunedMap, sequence)
    # print(internalNodes)
    # print(guestData)

    leafList = []
    marginal = open("/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_marginalAncestralStates/marginal_step_" + str(step), "w")
    leaves = open("/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_leaves/leaves_step_" + str(step) + ".fasta", "w")
    distance = cumulative(guest)
    distances = open("/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_distances/distances_step_" + str(step), "w")
    guestFile = open("/Users/williamlin/Desktop/IW/IW/phyloSim-master/simulated_guestTrees/guestTrees_step_" + str(step), "w")
    for node in internalNodes:
        search = guest&node
        distances.write(node + " " + str(distance[node]) + "\n")
        if search.is_leaf() == False:
            marginal.write(node + " " + internalNodes[node] + "\n")
        else:
            leafList.append(node)
            leaves.write(">" + node + "\n" + internalNodes[node] + "\n")
    guestFile.write(lucky.write(format=1))
    marginal.close()
    leaves.close()
    distances.close()
    guestFile.close()

    # print(guestData)
    aligns = {}
    # print(aligns)
    for guestNode in guestData:
        if guestNode in leafList:
            start = guestData[guestNode][1]
            end = guestData[guestNode][2]
            # print(start)
            # print(end)
            if (start,end) not in aligns:
                aligns[(start, end)] = {guestNode: guestData[guestNode][0]}
            else:
                if guestNode not in aligns[(start,end)]:
                    aligns[(start,end)][guestNode] = guestData[guestNode][0]

    increment = 0
    # print(aligns)
    # print(leafList)
    for column in aligns:
        if len(aligns[column]) > 8:
            columnFile = open("/Users/williamlin/Desktop/IW/IW/phyloSim-master/data/column" + str(increment) + "_" + str(step) + ".fasta",  "w")
            increment += 1
            for entry in aligns[column]:
                columnFile.write(">" + entry + "\n" + aligns[column][entry] + "\n")
            columnFile.close()

            



