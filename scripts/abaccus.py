#!/bin/python

'''Abaccus will incorporate taxonomic information in order to infer the
minimal number of losses that would explain a given taxonomic distribution
assuming only vertical inheritance.'''

from ete3 import Tree
import sys
import re
import os.path
import argparse
# from os.path import dirname (never used)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', default=False, required=True, help="Folder of a phylomizer output. It is expected to include a phylogenetic tree in newick format used as the input, and a rank file with the score for the different models. It will select the best available newick tree based on the score.")
	parser.add_argument('-m', '--main', help="This is the name of the central sequence used for building the tree. By default, the script will assume that the name of the sequence is included in the name of the file containing the tree and will search it using a regular expression match.")
	parser.add_argument('-o', '--output', default="abaccus_output.txt", help="Name of the output file")
	parser.add_argument('-t', '--taxonomy', default="./taxonomy.csv", help="Name of the file containing the taxonomic information of each species in the dataset")
	parser.add_argument('-J', '--jumps', default=2, type=int, help="Minimum number of unshared taxa between one branch of the tree and its sister branch, to be considered suspicious of beeing en event")
	parser.add_argument('-L', '--losses', default=3, type=int, help="Number of losses needed to be considered a potential event.")
	parser.add_argument('-R', '--rounds', default=6, type=int, help="It sets the number 'r' of times the algorithm will explore parent nodes. If they iterate r times without finding anything suspicious, the program will stop. This parameter aims to prevent false positives by limiting the explored space to the closest relatives to the central sequence.")
	parser.add_argument('-v', '--verbose', action='store_true', help="Prints the result onscreen instead of redirecting it to a file")
	parser.add_argument('-b', '--branch_support', default=0.9, type=float)
	# parser.add_argument('-X', '--exceptions', default=2, type=int, help="Number of acceptable consecutive breaks of monophily. The program will sample the tree when it arrives to this number. This number must be necessarily equal or lower than 'exceptions'.")
	# parser.add_argument('-S', '--strikes', default=3, type=int, help="Number of acceptable non-consecutive breaks of monophily. The program will sample the tree when it arrives to this number. This number must be necessarily equal or higher than 'exceptions'.")

args = parser.parse_args()

if args.input[-1] == "/":
	args.input = args.input[:-1]

# if args.strikes < args.exceptions:
# 	print("Error: Strikes must be equal or lower than Exceptions. Remember that	every Exception is an strike. but not every strike is an exception!")
# 	sys.exit()

if args.branch_support > 1.0:
	print("Branch support is a fraction of 1, so you cannot give a value above that. Your value was: ") + str(args.branch)
	sys.exit()

if os.path.isdir(args.input):
	for element in os.listdir(args.input):
		if element.find(".tree.phyml.rank.") > -1:
			rank = args.input + "/" + element
			break

	rankfile = open(rank)
	model = str(rankfile.readline()).split()[0]
	tree = rank[rank.rfind("/")+1:rank.find(".tree.phyml.rank")]
	attempt1 = args.input + "/" + tree + ".tree.phyml.ml." + model + ".nw"
	attempt2 = args.input + "/" + tree + ".tree.phyml.nj." + model + ".nw"
	chosen_model_tree = ''

	if os.path.isfile(attempt1):
		chosen_model_tree = attempt1
	else:
		if os.path.isfile(attempt2):
			chosen_model_tree = attempt2
		else:
			print("Wrong named or non-existent tree")
			sys.exit()

	phylotree = Tree(chosen_model_tree)
else:
	phylotree = Tree(args.input)

for node in phylotree.get_descendants():
	if node.support < args.branch_support:
		node.delete()

if args.main is not None:
	central_seq = args.main
else:
	if os.path.isdir(args.input):
		tag = str(args.input[args.input.rfind("."):])
		central_seq = re.findall("[A-Z0-9]{4,10}_[A-Z0-9]{3,5}" + tag, args.input)[-1][:-6]
	else:
		for seq in phylotree.get_leaves():
			name = str(seq)[3:str(seq).find("|")]
			if args.input.find(name) > -1:
				central_seq = re.findall("[A-Z0-9]{4,10}_[A-Z0-9]{3,5}", str(seq))[-1]
				break
		else:
			raise IOError("I don't know which one is the main sequence used for the tree. Maybe you can provide the name manually with the -m or --main flag")
	print('the central sequence found is: ' + str(central_seq))

central_sp = central_seq[central_seq.find("_")+1:]
taxofile = args.taxonomy
if not args.verbose:
	outputfile = open(args.output, "a")

########################################################################
'''Taxonomist will create a dictionary relating each mnemonic with the
complete taxonomy of this organism'''

taxo_dict = {}


def taxonomist(taxofile, sep=";"):
	for line in open(taxofile):
		if line[0] == "#":
			continue
		# strip newline
		line = line.strip()
		# skip first header line
		taxa = line.split(sep)
		# requires second column to be the mnemo
		taxo_dict[taxa[1]] = taxa[2:]
	return taxo_dict


########################################################################
'''Paperbag will create a dictionary relating each taxonomic term with
the species that are included in'''

paperbag_dict = {}


def paperbag(taxofile, sep=";"):
	# init empty set
	termset = set()
	# add most basal key (whichever is not EUKA!)
	paperbag_dict["biosphere"] = set()
	for line in open(taxofile):
		if line[0] == "#":
			continue
		else:
			line = line.strip()
			cast = line.split(sep)
			# assume third col is the first useful
			for term in cast[2:]:
				termset.add(term)
	for term in termset:
		# for each term in taxofile add a key and init empty set
		paperbag_dict[term] = set()
	for line in open(taxofile):
		if line[0] == "#":
			continue
		else:
			line = line.strip()
			elements = line.split(sep)
			# for each element add the corresponding mnemos
			for element in elements[2:]:
				paperbag_dict[element].add(elements[1])
				paperbag_dict["biosphere"].add(elements[1])
	return paperbag_dict


########################################################################
'''Dinasty will search the common ancestor of all the mnemonics in a list'''


def dinasty(spp_list):
	spp_tuple = tuple(spp_list)
	commonancestor = ''
	checkset = set()
	switch = False
	set_len = {len(taxo_dict[k]) for k in taxo_dict.keys()}
	if len(set_len) == 1:
		level = tuple(set_len)[0]
	else:
		print("Warning: Different taxonomy sizes")
	for i in range(level):
		if switch:
			continue
		else:
			for spp in spp_tuple:
				if spp in taxo_dict.keys():
					checkset.add(taxo_dict[spp][int(i)])
				if spp not in taxo_dict.keys():
					commonancestor = ["biosphere", level+1]
					switch = True
					break
		if len(checkset) == 1 and not switch:
			commonancestor = [taxo_dict[spp_list[0]][int(i)], int(i)]
			break
		else:
			checkset = set()
			continue
	if commonancestor == '':
		commonancestor = ['Eukaryota', level]
	return commonancestor


########################################################################
'''Get_mnemonics will use a phylotree or subset of it as input and will return
a list containing all the mnemonics included in it'''


def get_mnemonics(branch):
	mnemoset = {str(leaf)[str(leaf).rfind("_")+1:] for leaf in branch.get_leaves()}
	mnemolist = [element for element in mnemoset]
	return mnemolist


########################################################################
'''Orthogroup will look at the sequence that is referred in the name of
the file. Then it will start climbing and checking the common ancestor
between the first sequence and the sister branches. When the common ancestor
of any of this sister branches is closer to the main sequence than the
previous sister branch (a.k.a broken monophily) it will stop and will
select the monophiletic subtree for further analysis.
Be noted that the name of the sequences contains the mnemonic of each
species, and thus must be extracted from there'''


def orthogroup(phylotree):
	main_leaf = ''
	for leaf in phylotree.get_leaves():
		if str(leaf).find(central_seq) > -1:
			main_leaf = leaf
			break
	# main_leaf = phylotree&central_seq
	tree_root = main_leaf.get_farthest_node()[0]
	phylotree.set_outgroup(tree_root)
	output = main_leaf
	start_point = main_leaf
	switch = False
	leap = 0
	rounds = 0
	while rounds <= args.rounds and not switch:
		sis_branch = start_point.get_sisters()
		next_branch = start_point.get_common_ancestor(sis_branch).get_sisters()
		start_mnemo = get_mnemonics(start_point)
		# All this loops will stablish the phylogenetic relationships between the branches analyzed
		sis_mnemo = [spp for element in sis_branch for spp in get_mnemonics(element)]

		next_mnemo = [spp for element in next_branch for spp in get_mnemonics(element)]

		# This step creates sets that combines start_point + sis_branch (combined) and sis_branch + next_branch (para)
		# This loop is needed so start_mnemo is independent from combined mnemo
		combined_mnemo = [spp for element in start_point for spp in get_mnemonics(element)] + [el for el in sis_mnemo]

		# This loop is needed so para_mnemo is independent from next_mnemo
		para_mnemo = [spp for element in next_branch for spp in get_mnemonics(element)] + [el for el in sis_mnemo]

		start_dinasty, combined_dinasty, para_dinasty, next_dinasty = dinasty(start_mnemo), dinasty(combined_mnemo), dinasty(para_mnemo), dinasty(next_mnemo)

		# if at least one in both start+sis and next is proka, check if all are proka. If not, go to next branch. If yes, get common ancestor as result
		if combined_dinasty[0] == "biosphere" and next_dinasty[0] == "biosphere" and (combined_dinasty[1] - start_dinasty[1]) >= args.jumps:  # This means that at least one sequence in both sister_branch and next_branch is prokaryotic
			bacteria_purity = True
			for element in sis_mnemo:
				if element in taxo_dict.keys():
					bacteria_purity = False
			for element in next_mnemo:
				if element in taxo_dict.keys():
					bacteria_purity = False
			if not bacteria_purity:
				start_point = start_point.get_common_ancestor(sis_branch)
				rounds = rounds + 1
				continue
			if bacteria_purity:
				output = start_point.get_common_ancestor(next_branch)
				switch = True
				leap = combined_dinasty[1] - start_dinasty[1]
				break

		elif combined_dinasty[0] != "biosphere" and next_dinasty[0] != "biosphere" and (combined_dinasty[1] - start_dinasty[1]) >= args.jumps:  # This means that no prokaryotic sequence have been found so far neither in sis_branch nor in next_branch, but the taxonomic jump parameter has been met. It indicates a possible case of intereukaryotic HGT
			# if not all the species contained in the start dinasty are into the para (at least one loss) then exit and get the common ancestor as result
			if not paperbag_dict[start_dinasty[0]].issubset(paperbag_dict[para_dinasty[0]]):
				output = start_point.get_common_ancestor(next_branch)
				switch = True
				leap = combined_dinasty[1] - start_dinasty[1]
				break
			else:
				start_point = start_point.get_common_ancestor(sis_branch)
				rounds = rounds + 1
		else:
			start_point = start_point.get_common_ancestor(sis_branch)
			rounds = rounds + 1
	return output, start_point, leap

########################################################################


'''Abaccus will check the tree extracted from orthogroup and will count
the amount of gene losses that are needed to explain the given taxonomic
distribution'''


def abaccus(suspect_tree, start_point):
	# takes the orthogroup output and analyse it by:
	observed_mnemo_list = set()
	event_mnemo_list = []
	wanted = []
	# observed mnemo has the leaves in the starting node
	for leaf in start_point.get_leaves():
		observed_mnemo_list.add(str(str(leaf)[str(leaf).rfind("_")+1:]))
	# event mnemo has all the leaf in suspect tree
	for leaf in suspect_tree.get_leaves():
		event_mnemo_list.append(str(str(leaf)[str(leaf).rfind("_")+1:]))
	# find common ancestor between nodes in starting point
	commonancestor = dinasty(event_mnemo_list)
	losses = 0
	# loop from 0 to Number of commonancestor rank (last level excluded)
	for i in range(commonancestor[1]):
		if i > len(taxo_dict[central_sp]) - 1:
			losses = losses + 1
			wanted.append("There is at least one loss in the base of eukarya")
			break
		# elif the species in the ith level are a subset of the species in the starting node
		elif set(paperbag_dict[taxo_dict[central_sp][i]]).issubset(observed_mnemo_list):
			continue
		# elif the length of the difference between ith level and the previous is 0
		elif len(set(paperbag_dict[taxo_dict[central_sp][i]]).difference(set(paperbag_dict[taxo_dict[central_sp][i-1]]))) == 0:
			continue
		else:
			losses = losses + 1
			wanted.append("There is at least one loss in " + taxo_dict[central_sp][i])
			continue
	return losses, wanted

########################################################################

taxonomist(taxofile)
paperbag(taxofile)
orthotree = orthogroup(phylotree)
jump = str(orthotree[2])
cutoff = abaccus(orthotree[0], orthotree[1])

if cutoff[0] >= args.losses:
	if not args.verbose:
		outputfile.write("\n###" + central_seq + ' - ' + taxo_dict[central_sp][0] + "###\n")
		outputfile.write("Minimal number of losses: " + str(cutoff[0]) + "\n")
		outputfile.write("J == " + jump + " ; L == " + str(cutoff[0]) + ". Cutoff values are J =< " + str(args.jumps) + " and L =< " + str(args.losses) + "\n")
		for element in cutoff[1]:
			outputfile.write(element)
			outputfile.write("\n")
		outputfile.write(orthotree[0].get_ascii(show_internal=False))
		outputfile.write("\n\n\n")
		outputfile.close()
	if args.verbose:
		print("###" + central_seq + "###")
		print("Minimal number of losses: " + str(cutoff[0]))
		print("J == " + jump + " and L == " + str(cutoff[0]) + ". Cutoff values are J => " + str(args.jumps) + " and L => " + str(args.losses) + "\n")
		for element in cutoff[1]:
			print(element)
		print(orthotree[0])
		print("\n")
else:
	print('Minimal number of losses: ' + str(cutoff[0]))
	print("J == " + jump + " and L == " + str(cutoff[0]) + ". Cutoff values are J => " + str(args.jumps) + " and L => " + str(args.losses))
	print('No events detected, empty output file')
