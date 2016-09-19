#!/bin/python

'''Abaccus will incorporate taxonomic information in order to infer the 
minimal number of loses that would explain a given taxonomic distribution
assuming only vertical inheritance.'''

from ete2 import Tree
import sys, re, os.path
import argparse
from os.path import dirname

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default=False, required=True, help="Folder of a phylomizer output. It is expected to include a phylogenetic tree in newick format used as the input, and a rank file with the score for the different models. It will select the best available newick tree based on the score.")
    parser.add_argument('-m', '--main', help="This is the name of the central sequence used for building the tree. By default, the script will assume that the name of the sequence is included in the name of the file containing the tree and will search it using a regular expression match.")
    parser.add_argument('-o', '--output', default="abaccus_output.txt", help="Name of the output file")
    parser.add_argument('-t', '--taxonomy', default="./taxonomy.csv", help="Name of the file containing the taxonomic information of each species in the dataset")
    parser.add_argument('-J', '--jumps', default=3, type=int, help="Minimum number of unshared taxa between one branch of the tree and its sister branch, to be considered suspicious of beeing en event")
    parser.add_argument('-L', '--loses', default=4, type=int, help="Number of loses needed to be considered a potential event.")
    parser.add_argument('-R', '--rounds', default=6, type=int, help="It sets the number 'r' of times the algorythm will explore parent nodes. If they iterate r times without finding anything suspicious, the program will stop. This parameter aims to prevent false positives by limiting the explored space to the closest relatives to the central sequence.")
    parser.add_argument('-X', '--exceptions', default=2, type=int, help="Number of acceptable consecutive breaks of monophily. The program will sample the tree when it arrives to this number. This number must be necessarily equal or lower than 'exceptions'.")
    parser.add_argument('-S', '--strikes', default=3, type=int, help="Number of acceptable non-consecutive breaks of monophily. The program will sample the tree when it arrives to this number. This number must be necessarily equal or higher than 'exceptions'.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Prints the result onscreen instead of redirecting it to a file")
    parser.add_argument('-b', '--branch_support', default=0.9, type=float)
        
args = parser.parse_args()

if args.input[-1] == "/":
	args.input = args.input[:-1]

if args.strikes < args.exceptions:
	print "Error: Strikes must be equal or lower than Exceptions. Remember that every Exception is an strike. but not every strike is an exception!"
	sys.exit()
	
if args.branch_support > 1.0:
	print "Branch support is a fraction of 1, so you cannot give a value above that. Your value was: "+str(args.branch)
	sys.exit()
  


if os.path.isdir(args.input) == True:
	for element in os.listdir(args.input):
		if element.find(".tree.phyml.rank.nj") > -1:
			rank = args.input + "/" + element
			break 			
	rankfile = open(rank)

	model = str(rankfile.readline()).split()[0]

	tree = rank[rank.rfind("/")+1:rank.find(".tree.phyml.rank")]
	attempt1 = args.input + "/" + tree + ".tree.phyml.ml." + model + ".nw"
	attempt2 = args.input + "/" + tree + ".tree.phyml.nj." + model + ".nw"

	chosen_model_tree = ''

	if os.path.isfile(attempt1) == True:
		chosen_model_tree = attempt1
	else:
		if os.path.isfile(attempt2) == True:
			chosen_model_tree = attempt2
		else:
			print "Wrong named or non-existent tree"
			sys.exit()
		
	phylotree = Tree(chosen_model_tree)
else:
	phylotree = Tree(args.input)

for node in phylotree.get_descendants():
	if node.support < args.branch_support:
		node.delete()
tag = str(args.input[args.input.rfind("."):])


if args.main != None:
	central_seq = args.main
else:
	if os.path.isdir(args.input) == True:
		central_seq = re.findall("[A-Z0-9]{4,10}_[A-Z0-9]{3,5}"+tag, args.input)[-1][:-6]
	else:
		for seq in phylotree.get_leaves():
			name = str(seq)[3:str(seq).find("|")]
			if args.input.find(name) > -1:
				central_seq = re.findall("[A-Z0-9]{4,10}_[A-Z0-9]{3,5}", str(seq))[-1]
				break
		else:
			raise IOError ("I don't know which one is the main sequence used for the tree. Maybe you can provide the name manually with the -m or --main flag")


central_sp = central_seq[central_seq.find("_")+1:]
taxofile = args.taxonomy
if args.verbose == False:
	outputfile = open(args.output, "a")

########################################################################
'''Taxonomist will create a dictionary relating each mnemonic with the
complete taxonomy of this organism'''

taxo_dict = {}

def taxonomist (taxofile):
	for line in open(taxofile):
		taxa = re.split('[,;]', line)
		taxo_dict[taxa[1]] = taxa[2:10]
	return taxo_dict

########################################################################
'''Paperbag will create a dictionary relating each taxonomic term with 
the species that are included in'''

paperbag_dict = {}

def paperbag (taxofile):
	termset = set()
	paperbag_dict["biosphere"] = set()
	paperbag_dict["biosphere"].add("BACTERIA")
	paperbag_dict["biosphere"].add("ARCHAEA")
	paperbag_dict["eukaryota"] = set()
	for line in open(taxofile):
		if line[0] == "#": continue
		else:
			cast = re.split('[,;]', line)
			for term in cast[2:10]:
				termset.add(term)
	for term in termset:
		paperbag_dict[term] = set()
	for line in open(taxofile):
		if line[0] == "#": continue
		else:
			elements = re.split('[,;]', line)
			for element in elements[2:10]:
				paperbag_dict[element].add(elements[1])
				paperbag_dict["biosphere"].add(elements[1])
				paperbag_dict["eukaryota"].add(elements[1])
	return paperbag_dict

########################################################################
'''euka_exclusive will look at a phylogenetic tree and check if all 
branches are eukaryotic, and some related information'''

def euka_exclusive (problem):
	cast = set()
	for leaf in problemtree.get_leaves():
		cast.add(str(leaf)[str(leaf).rfind("_")+1:])
	return cast.issubset(paperbag_dict["eukaryota"]), cast, cast.intersection(paperbag_dict["eukaryota"])

########################################################################
'''Dinasty will search the common ancestor of all the mnemonics in a list'''

def dinasty (spp_list):
	spp_tuple = tuple(spp_list)
	commonancestor = ''
	checkset = set()
	switch = False
	for i in range(len(taxo_dict["ARATH"])):
		if switch == True: continue
		else:
			for spp in spp_tuple:
				if spp in taxo_dict.keys():
					checkset.add(taxo_dict[spp][int(i)])
				if spp not in taxo_dict.keys():
					commonancestor = ["biosphere", 9]
					switch = True
					break
		if len(checkset) == 1 and switch == False:
			commonancestor = [taxo_dict[spp_list[0]][int(i)], int(i)]
			break
		else:
			checkset = set() 
			continue
	if commonancestor == '':
		commonancestor = ["eukaryota", 8]
	return commonancestor

########################################################################
'''Get_mnemonics will use a phylotree or subset of it as inpunt and will return
a list containing all the mnemonics included in it'''

def get_mnemonics(branch):
	mnemoset = set()
	mnemolist = []
	for leaf in branch.get_leaves():
		mnemoset.add(str(leaf)[str(leaf).rfind("_")+1:])
	for element in mnemoset:
		mnemolist.append(element)
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

def orthogroup (phylotree):
	main_leaf = ''
	for leaf in phylotree.get_leaves():
		if str(leaf).find(central_seq) > -1:
			main_leaf = leaf
			break
	distance_to_main = 0.0
	tree_root = ''			
	for leaf in phylotree.iter_leaves(): #With this, the script will look at the farthest node to the central sequence and will select it as the root of the tree
		if leaf.name != main_leaf:
			if distance_to_main < main_leaf.get_distance(leaf):
				distance_to_main = main_leaf.get_distance(leaf)
				tree_root = leaf
	phylotree.set_outgroup(tree_root) 		
	output = main_leaf
	prev_level = taxo_dict[central_sp][0]	
	curr_level = taxo_dict[central_sp][0]
	start_point = main_leaf
	switch = False
	rounds = 0
	while rounds <= args.rounds and switch == False:
		sis_mnemo, next_mnemo, combined_mnemo, para_mnemo, start_mnemo = [], [], [], [], []
		sis_branch = start_point.get_sisters()	
		next_branch = start_point.get_common_ancestor(sis_branch).get_sisters()
		start_mnemo = get_mnemonics(start_point)
		#All this loops will stablish the phylogenetic relationships between the branches analyzed
		for element in sis_branch:
			for spp in get_mnemonics(element):
				sis_mnemo.append(spp)
		
		for element in next_branch:
			for spp in get_mnemonics(element):
				next_mnemo.append(spp)
		
		for element in start_point: #This loop is needed so start_mnemo is independent from combined mnemo
			for spp in get_mnemonics(element):
				combined_mnemo.append(spp)
		
		for element in next_branch: #This loop is needed so para_mnemo is independent from next_mnemo
			for spp in get_mnemonics(element):
				para_mnemo.append(spp)
		
		
		for element in sis_mnemo: #This step creates sets that combines start_point + sis_branch (combined) and sis_branch + next_branch (para)
			combined_mnemo.append(element)
			para_mnemo.append(element)
			
		start_dinasty, combined_dinasty, para_dinasty, next_dinasty, sis_dinasty = dinasty(start_mnemo), dinasty(combined_mnemo), dinasty(para_mnemo), dinasty(next_mnemo), dinasty(sis_mnemo)

		if combined_dinasty[0] == "biosphere" and next_dinasty[0] == "biosphere" and (combined_dinasty[1] - start_dinasty[1]) >= args.jumps: #This means that at least one sequence in both sister_branch and next_branch is prokaryotic
			bacteria_purity = True
			for element in sis_mnemo:
				if element in taxo_dict.keys():
					bacteria_purity = False
			for element in next_mnemo:
				if element in taxo_dict.keys():
					bacteria_purity = False
			if bacteria_purity == False:
				start_point = start_point.get_common_ancestor(sis_branch)
				rounds = rounds + 1
				continue
			if bacteria_purity == True:
				output = start_point.get_common_ancestor(next_branch)
				switch = True
				break
			
		elif combined_dinasty[0] != "biosphere" and next_dinasty[0] != "biosphere" and (combined_dinasty[1] - start_dinasty[1]) >= args.jumps: #This means that no prokaryotic sequence have been found so far neither in sis_branch nor in next_branch, but the taxonomic jump parameter has been met. It indicates a possible case of intereukaryotic HGT
			if paperbag_dict[start_dinasty[0]].issubset(paperbag_dict[para_dinasty[0]]) == False:
				output = start_point.get_common_ancestor(next_branch)
				switch = True
				break
			else:
				start_point = start_point.get_common_ancestor(sis_branch)
				rounds = rounds + 1
		else:
			start_point = start_point.get_common_ancestor(sis_branch)
			rounds = rounds + 1
	return output, start_point, start_dinasty[0], sis_dinasty[0]

########################################################################
'''Abaccus will check the tree extracted from orthogroup and will count
the amount of gene loses that are needed to explain the given taxonomic 
distribution'''

def abaccus(suspect_tree, start_point, acceptor, donnor):
	observed_mnemo_list = set()
	event_mnemo_list = []
	wanted = []
	for leaf in start_point.get_leaves():
		observed_mnemo_list.add(str(str(leaf)[str(leaf).rfind("_")+1:]))
	for leaf in suspect_tree.get_leaves():
		event_mnemo_list.append(str(str(leaf)[str(leaf).rfind("_")+1:]))
	
	commonancestor = dinasty(event_mnemo_list)
	loses = 0
	for i in range(commonancestor[1]):
		if i > len(taxo_dict[central_sp]) - 1:
			loses = loses + 1
			wanted.append("There is at least one loss in the base of eukarya")
			break
		elif set(paperbag_dict[taxo_dict[central_sp][i]]).issubset(observed_mnemo_list):
			continue
		elif len(set(paperbag_dict[taxo_dict[central_sp][i]]).difference(set(paperbag_dict[taxo_dict[central_sp][i-1]]))) == 0:
			continue
		else:
			loses = loses +1
			wanted.append("There is at least one loss in " + taxo_dict[central_sp][i])
			continue
	return loses, wanted
				
########################################################################

taxonomist(taxofile)
paperbag(taxofile)
orthotree = orthogroup(phylotree)

cutoff = abaccus(orthotree[0], orthotree[1], orthotree[2], orthotree[3])

if cutoff[0] >= args.loses:
	if args.verbose == False:
		outputfile.write("\n###" + central_seq + "###\n")
		outputfile.write("#Minimal number of loses: " + str(cutoff[0]) + "\n")
		for element in cutoff[1]:
			outputfile.write(element)
			outputfile.write("\n")
		outputfile.write(orthotree[0].get_ascii(show_internal=False))
		outputfile.write("\n\n\n")
		outputfile.close()
	if args.verbose == True:
		print "###"+central_seq+"###"
		print("#Minimal number of loses: " + str(cutoff[0]))
		for element in cutoff[1]:
			print(element)
		print(orthotree[0])
		print ("\n")
