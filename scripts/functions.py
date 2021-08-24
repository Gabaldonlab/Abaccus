from ete3 import Tree, PhyloTree
import sys
import re
import os.path


def taxonomist(taxofile, phylo=False, sep=";", proka=None):

    '''Taxonomist will create a dictionary relating each mnemonic with the
    complete taxonomy of this organism'''

    taxo_dict = {}

    if phylo:
        num = 1
        # name intenal nodes with arbitrary nodes
        for node in taxofile.iter_descendants():
            if node.name == "":
                node.name = "internal_" + str(num)
                num += 1
        for node in taxofile.iter_leaves():
            value = [nod.name for nod in node.get_ancestors()][:-1]
            # the user decides if excluding proka entries
            if proka is not None:
                # only non proka entries will be keys
                if node.name not in proka:
                    taxo_dict[node.name] = value
            else:
                taxo_dict[node.name] = value
    else:
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


def paperbag(taxofile, phylo=False, sep=";", proka=None):

    '''Paperbag will create a dictionary relating each taxonomic term with
    the species that are included in'''

    paperbag_dict = {}

    if phylo:
        for node in taxofile.iter_descendants():
            # if node.name in taxo_dict.keys():
            mnemos = set(node.get_leaf_names())
            if proka is not None:
                # only add set without proka entries
                if len(mnemos.intersection(set(proka)))==0:
                    paperbag_dict[node.name] = mnemos
            else:
                paperbag_dict[node.name] = mnemos
        # add all to root and refer to root as biosphere
        if proka is not None:
            paperbag_dict['biosphere'] = set(taxofile.get_leaf_names()) - set(proka)
        else:
            paperbag_dict['biosphere'] = set(taxofile.get_leaf_names())
    else:
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


def dinasty(spp_list, taxo_dict):

    '''Dinasty will search the common ancestor of all the mnemonics in a list'''

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


def dinasty_tree(spp_list, sp_tree, taxo_dict, spe2age):
    # if (sp_tree.get_common_ancestor(spp_list) == sp_tree.get_tree_root()
    #     and len(spp_list) > 1):
    if len(set(spp_list) - set(taxo_dict.keys())) > 0:
        # biosphere (not proper if all species are euka or proka) if
        commonancestor = ["biosphere", max(spe2age.values())]
    else:
        if len(spp_list) > 1:
            anc = sp_tree.get_common_ancestor(spp_list)
        else:
            anc = sp_tree.get_leaves_by_name(spp_list[0])[0]
        ages = [spe2age[k] for k in anc.get_leaf_names()]
        val = max(ages)  # - min(ages)
        commonancestor = [anc.name, val]
    return commonancestor


########################################################################


def get_mnemonics(branch):
    '''Get_mnemonics will use a phylotree or subset of it as input and will return
    a list containing all the mnemonics included in it'''
    mnemoset = {str(leaf)[str(leaf).rfind("_")+1:] for leaf in branch.get_leaves()}
    mnemolist = [element for element in mnemoset]
    return mnemolist


########################################################################


def load_species_name(node):
    return node.split("_")[1]

def root_tree(phylotree, main_leaf, spe2age=None, midpoint=True):
    if midpoint:
        tree_root = main_leaf.get_farthest_node()[0]
        phylotree.set_outgroup(tree_root)
    if not midpoint and spe2age is not None:
        phylotree = PhyloTree(phylotree.write(), sp_naming_function=load_species_name)
        phylotree.set_outgroup(phylotree.get_farthest_oldest_leaf(spe2age))

    return phylotree


def orthogroup(phylotree, taxo_dict, paperbag_dict, main_leaf, rounds_def, jumps):
    '''Orthogroup will look at the sequence that is referred in the name of
    the file. Then it will start climbing and checking the common ancestor
    between the first sequence and the sister branches. When the common ancestor
    of any of this sister branches is closer to the main sequence than the
    previous sister branch (a.k.a broken monophily) it will stop and will
    select the monophiletic subtree for further analysis.
    Be noted that the name of the sequences contains the mnemonic of each
    species, and thus must be extracted from there'''

    output = main_leaf
    start_point = main_leaf
    switch = False
    leap = 0
    rounds = 0
    while rounds <= rounds_def and not switch:
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

        start_dinasty, combined_dinasty, para_dinasty, next_dinasty = dinasty(start_mnemo, taxo_dict), dinasty(combined_mnemo, taxo_dict), dinasty(para_mnemo, taxo_dict), dinasty(next_mnemo, taxo_dict)

        # if at least one in both start+sis and next is proka, check if all are proka. If not, go to next branch. If yes, get common ancestor as result
        if combined_dinasty[0] == "biosphere" and next_dinasty[0] == "biosphere" and (combined_dinasty[1] - start_dinasty[1]) >= jumps:  # This means that at least one sequence in both sister_branch and next_branch is prokaryotic
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

        elif combined_dinasty[0] != "biosphere" and next_dinasty[0] != "biosphere" and (combined_dinasty[1] - start_dinasty[1]) >= jumps:  # This means that no prokaryotic sequence have been found so far neither in sis_branch nor in next_branch, but the taxonomic jump parameter has been met. It indicates a possible case of intereukaryotic HGT
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


def orthogroup_tree(phylotree, taxo_dict, paperbag_dict, main_leaf, rounds_def, jumps, ref_tree, sp2age):
    output = main_leaf
    start_point = main_leaf
    switch = False
    leap = 0
    rounds = 0
    while rounds <= rounds_def and not switch:
        sis_branch = start_point.get_sisters()
        next_branch = start_point.get_common_ancestor(sis_branch).get_sisters()
        start_mnemo = get_mnemonics(start_point)
        # All this loops will stablish the phylogenetic relationships between the branches analyzed
        sis_mnemo = [spp for element in sis_branch for spp in get_mnemonics(element)]

        next_mnemo = [spp for element in next_branch for spp in get_mnemonics(element)]

        # This step creates sets that combines start_point + sis_branch (combined) and sis_branch + next_branch (para)
        # This loop is needed so start_mnemo is independent from combined mnemo
        combined_mnemo = [
            spp for element in start_point for spp in get_mnemonics(element)
        ] + [el for el in sis_mnemo]

        # This loop is needed so para_mnemo is independent from next_mnemo
        para_mnemo = [
            spp for element in next_branch for spp in get_mnemonics(element)
        ] + [el for el in sis_mnemo]

        start_dinasty, combined_dinasty, para_dinasty, next_dinasty = (
            dinasty_tree(start_mnemo, ref_tree, taxo_dict, sp2age),
            dinasty_tree(combined_mnemo, ref_tree, taxo_dict, sp2age),
            dinasty_tree(para_mnemo, ref_tree, taxo_dict, sp2age),
            dinasty_tree(next_mnemo, ref_tree, taxo_dict, sp2age),
        )

        # if at least one in both start+sis and next is proka, check if all are proka. If not, go to next branch. If yes, get common ancestor as result
        if (
            combined_dinasty[0] == "biosphere"
            and next_dinasty[0] == "biosphere"
            and (combined_dinasty[1] - start_dinasty[1]) >= jumps
        ):  # This means that at least one sequence in both sister_branch and next_branch is prokaryotic
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

        elif (
            combined_dinasty[0] != "biosphere"
            and next_dinasty[0] != "biosphere"
            and (combined_dinasty[1] - start_dinasty[1]) >= jumps
        ):  # This means that no prokaryotic sequence have been found so far neither in sis_branch nor in next_branch, but the taxonomic jump parameter has been met. It indicates a possible case of intereukaryotic HGT
            # if not all the species contained in the start dinasty are into the para (at least one loss) then exit and get the common ancestor as result
            if not paperbag_dict[start_dinasty[0]].issubset(
                paperbag_dict[para_dinasty[0]]
            ):
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


def abaccus(suspect_tree, start_point, taxo_dict, paperbag_dict, central_sp):
    '''Abaccus will check the tree extracted from orthogroup and will count
    the amount of gene losses that are needed to explain the given taxonomic
    distribution'''
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
    commonancestor = dinasty(event_mnemo_list, taxo_dict)
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


def abaccus_tree(suspect_tree, start_point, taxo_dict, paperbag_dict, ref_tree, sp2age, central_sp):
    pass
    # this still does not make sens with a phylogeny
    # observed_mnemo_list = set()
    # event_mnemo_list = []
    # wanted = []
    # # observed mnemo has the leaves in the starting node
    # for leaf in start_point.get_leaves():
    #     observed_mnemo_list.add(str(str(leaf)[str(leaf).rfind("_") + 1 :]))
    # # event mnemo has all the leaf in suspect tree
    # for leaf in suspect_tree.get_leaves():
    #     event_mnemo_list.append(str(str(leaf)[str(leaf).rfind("_") + 1 :]))
    # # find common ancestor between nodes in starting point
    # commonancestor = dinasty_tree(event_mnemo_list, ref_tree, taxo_dict, sp2age)
    # losses = 0
    # # loop from 0 to Number of commonancestor rank (last level excluded)
    # for i in range(commonancestor[1]):
    #     if i > len(taxo_dict[central_sp]) - 1:
    #         losses = losses + 1
    #         wanted.append("There is at least one loss in the base of eukarya")
    #         break
    #         # elif the species in the ith level are a subset of the species in the starting node
    #     elif set(paperbag_dict[taxo_dict[central_sp][i]]).issubset(
    #         observed_mnemo_list
    #     ):
    #         continue
    #     # elif the length of the difference between ith level and the previous is 0
    #     elif (
    #         len(
    #             set(paperbag_dict[taxo_dict[central_sp][i]]).difference(
    #                 set(paperbag_dict[taxo_dict[central_sp][i - 1]])
    #             )
    #         )
    #         == 0
    #     ):
    #         continue
    #     else:
    #         losses = losses + 1
    #         wanted.append(
    #             "There is at least one loss in " + taxo_dict[central_sp][i]
    #         )
    #         continue
    #
    # return losses, wanted


#########################################################################


def get_sp2age(sptree, central_sp):
    main_leaf = sptree & central_sp
    # compute sp2age from species tree
    sp2age = {}
    for leaf in sptree.iter_leaves():
        mrca = main_leaf.get_common_ancestor(leaf)
        d = main_leaf.get_distance(mrca, topology_only=True) + 1
        sp2age[leaf.name] = int(d)
    return sp2age
