#!/bin/python

"""Abaccus will incorporate taxonomic information in order to infer the
minimal number of losses that would explain a given taxonomic distribution
assuming only vertical inheritance."""

import sys

sys.path.insert(1, "scripts")
from ete3 import Tree, TreeStyle
import re
import os
import argparse
import functions as fn

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default=False,
        required=True,
        help="a **phylogenetic tree** representing the evolution of a gene family. The **phylogenetic tree** can be provided in standard **newick** format, or as provided by the **phylomizer pipeline** (https://github.com/Gabaldonlab/phylomizer). In that case the input should be the folder of a phylomizer output. It is expected to include a phylogenetic tree in newick format used as the input, and a rank file with the score for the different models. It will select the best available newick tree based on the score.",
    )
    parser.add_argument(
        "-m",
        "--main",
        help="This is the name of the central sequence used for building the tree. By default, the script will assume that the name of the sequence is included in the name of the file containing the tree and will search it using a regular expression match.",
    )
    parser.add_argument(
        "-o", "--output", default="abaccus_output.txt", help="Name of the output file"
    )
    parser.add_argument(
        "-t",
        "--taxonomy",
        help="The true distribution of species either as: a **.csv** file containing a **taxonomy** (i.e data/taxonomy.csv, included in this folder). The **taxonomy** should include at least all species contained in the **phylogenetic tree**. Or a species tree where the true phyloegenetic relationships between all the species present in the gene family tree are present. This file must end with **.newick**,**.nw** or **.nwk**.",
    )
    parser.add_argument(
        "-n",
        "--naming",
        default="pdb",
        type=str,
        help="Naming scheme for species annotation in the gene tree. Value can be 'pdb' for phylome db style or 'uniref' for uniref style",
    )
    parser.add_argument(
        "-J",
        "--jumps",
        default=2,
        type=int,
        help="Minimum number of unshared taxa between one branch of the tree and its sister branch, to be considered suspicious of being en event",
    )
    parser.add_argument(
        "-L",
        "--losses",
        default=3,
        type=int,
        help="Number of losses needed to be considered a potential event.",
    )
    parser.add_argument(
        "-R",
        "--rounds",
        default=6,
        type=int,
        help="It sets the number 'r' of times the algorithm will explore parent nodes. If they iterate r times without finding anything suspicious, the program will stop. This parameter aims to prevent false positives by limiting the explored space to the closest relatives to the central sequence.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Prints the result onscreen instead of redirecting it to a file",
    )
    parser.add_argument(
        "--image",
        action="store_true",
        help="Produces the abaccus plot in working dir or in output directory if output is given (blue point is starting point, red point is end point)",
    )
    parser.add_argument("-b", "--branch_support", default=0.9, type=float)
    parser.add_argument(
        "-p",
        "--prokaryotes",
        help="Name of the file containing the mnemonic codes of prokaryotic entitites. Use only if using abaccus with a phylogeny.",
    )

    # parser.add_argument('-X', '--exceptions', default=2, type=int, help="Number of acceptable consecutive breaks of monophily. The program will sample the tree when it arrives to this number. This number must be necessarily equal or lower than 'exceptions'.")
    # parser.add_argument('-S', '--strikes', default=3, type=int, help="Number of acceptable non-consecutive breaks of monophily. The program will sample the tree when it arrives to this number. This number must be necessarily equal or higher than 'exceptions'.")


args = parser.parse_args()

if args.input[-1] == "/":
    args.input = args.input[:-1]

taxofile = args.taxonomy

if args.naming not in ["pdb", "uniref"] and not "":
    sys.exit("naming scheme must be either pdb or uniref")

naming_scheme = args.naming


if taxofile.endswith("csv"):
    phylogeny = False
elif taxofile.endswith(("nwk", "nw", "newick")):
    phylogeny = True
    species_tree = Tree(taxofile)
    if args.prokaryotes:
        with open(args.prokaryotes) as p:
            proka = [line.strip() for line in p.readlines()]
    else:
        proka = None
else:
    sys.exit(
        "Please end the taxonomy file with csv or the phylogeny with nw, nwk or newick in order for the program to run accordingly"
    )


# if args.strikes < args.exceptions:
#     print("Error: Strikes must be equal or lower than Exceptions. Remember that    every Exception is an strike. but not every strike is an exception!")
#     sys.exit()

if args.branch_support > 1.0:
    print(
        "Branch support is a fraction of 1, so you cannot give a value above that. Your value was: "
    ) + str(args.branch)
    sys.exit()

if os.path.isdir(args.input):
    for element in os.listdir(args.input):
        if element.find(".tree.phyml.rank.") > -1:
            rank = args.input + "/" + element
            break
    rankfile = open(rank)
    model = str(rankfile.readline()).split()[0]
    tree = rank[rank.rfind("/") + 1 : rank.find(".tree.phyml.rank")]
    attempt1 = args.input + "/" + tree + ".tree.phyml.ml." + model + ".nw"
    attempt2 = args.input + "/" + tree + ".tree.phyml.nj." + model + ".nw"
    chosen_model_tree = ""
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

phylotree = fn.load_tree(phylotree, naming_scheme=naming_scheme)

for node in phylotree.get_descendants():
    if node.support < args.branch_support:
        node.delete()

if args.main is not None:
    central_seq = args.main
else:
    if os.path.isdir(args.input):
        tag = str(args.input[args.input.rfind(".") :])
        central_seq = re.findall("[A-Z0-9]{4,10}_[A-Z0-9]{3,5}" + tag, args.input)[-1][
            :-6
        ]
    else:
        for seq in phylotree.get_leaves():
            name = str(seq)[3 : str(seq).find("|")]
            if args.input.find(name) > -1:
                central_seq = re.findall("[A-Z0-9]{4,10}_[A-Z0-9]{3,5}", str(seq))[-1]
                break
        else:
            raise IOError(
                "I don't know which one is the main sequence used for the tree. Maybe you can provide the name manually with the -m or --main flag"
            )
    print("the central sequence found is: " + str(central_seq))


main_leaf = ""
for leaf in phylotree.get_leaves():
    if str(leaf).find(central_seq) > -1:
        main_leaf = leaf
        break

central_sp = main_leaf.species


out_img = ""
if not args.verbose:
    outputfile = open(args.output, "a")
    outfile_noext = os.path.basename(args.output).split(".")[0]
    outputdir = os.path.dirname(args.output)
    out_img = outputdir + "/" + outfile_noext + "_abaccus.png"

# add option how to root it, now its midpoint
root_phylotree = fn.root_tree(phylotree, main_leaf)

if phylogeny:
    species_tree = fn.name_internal(species_tree)
    list_gene_sp = list(set([sp.species for sp in phylotree.get_leaves()]))
    species_tree.prune(list_gene_sp)
    taxo_dict = fn.taxonomist(species_tree, phylo=True, proka=proka)
    paperbag_dict = fn.paperbag(species_tree, phylo=True, proka=proka)
    spe2age = fn.get_sp2age(species_tree, central_sp)
    orthotree = fn.orthogroup_tree(
        root_phylotree,
        taxo_dict,
        paperbag_dict,
        main_leaf,
        args.rounds,
        args.jumps,
        species_tree,
        spe2age,
    )
    losses = fn.abaccus_tree(orthotree[0], species_tree, main_leaf)
    if losses >= args.losses:
        if not args.verbose:
            outputfile.write("\n###" + central_seq + " - " + central_sp + "###\n")
            outputfile.write("Minimal number of losses: " + str(losses) + "\n")
            outputfile.write(
                "J == "
                + orthotree[2]
                + " ; L == "
                + str(losses)
                + ". Cutoff values are J =< "
                + str(args.jumps)
                + " and L =< "
                + str(args.losses)
                + "\n"
            )
            # for element in cutoff[1]:
            #     outputfile.write(element)
            #     outputfile.write("\n")
            outputfile.write(orthotree[0].get_ascii(show_internal=False))
            outputfile.write("\n\n\n")
            outputfile.close()
        if args.verbose:
            print("###" + central_seq + "###")
            print("Minimal number of losses: " + str(losses))
            print(
                "J == "
                + str(orthotree[2])
                + " and L == "
                + str(losses)
                + ". Cutoff values are J => "
                + str(args.jumps)
                + " and L => "
                + str(args.losses)
                + "\n"
            )
            # for element in cutoff[1]:
            #     print(element)
            print(orthotree[0])
            print("\n")
        if args.image:
            if out_img == "":
                out_img = central_seq + "_abaccus.png"
            ts = fn.gene_viz(root_phylotree, orthotree, main_leaf)
            phylotree.render(out_img, tree_style=ts)
    else:
        print("Minimal number of losses: " + str(losses))
        print(
            "J == "
            + str(orthotree[2])
            + " and L == "
            + str(losses)
            + ". Cutoff values are J => "
            + str(args.jumps)
            + " and L => "
            + str(args.losses)
        )
        print("No events detected, empty output file")

else:
    taxo_dict = fn.taxonomist(taxofile)
    paperbag_dict = fn.paperbag(taxofile)
    orthotree = fn.orthogroup(
        root_phylotree, taxo_dict, paperbag_dict, main_leaf, args.rounds, args.jumps
    )
    jump = str(orthotree[2])
    cutoff = fn.abaccus(
        orthotree[0], orthotree[1], taxo_dict, paperbag_dict, central_sp
    )

    if cutoff[0] >= args.losses:
        if not args.verbose:
            outputfile.write(
                "\n###" + central_seq + " - " + taxo_dict[central_sp][0] + "###\n"
            )
            outputfile.write("Minimal number of losses: " + str(cutoff[0]) + "\n")
            outputfile.write(
                "J == "
                + jump
                + " ; L == "
                + str(cutoff[0])
                + ". Cutoff values are J =< "
                + str(args.jumps)
                + " and L =< "
                + str(args.losses)
                + "\n"
            )
            for element in cutoff[1]:
                outputfile.write(element)
                outputfile.write("\n")
            outputfile.write(orthotree[0].get_ascii(show_internal=False))
            outputfile.write("\n\n\n")
            outputfile.close()
        if args.verbose:
            print("###" + central_seq + "###")
            print("Minimal number of losses: " + str(cutoff[0]))
            print(
                "J == "
                + jump
                + " and L == "
                + str(cutoff[0])
                + ". Cutoff values are J => "
                + str(args.jumps)
                + " and L => "
                + str(args.losses)
                + "\n"
            )
            for element in cutoff[1]:
                print(element)
            print(orthotree[0])
            print("\n")
        if args.image:
            if out_img == "":
                out_img = central_seq + "_abaccus.png"
            ts = fn.gene_viz(root_phylotree, orthotree, main_leaf)
            phylotree.render(out_img, tree_style=ts)
    else:
        print("Minimal number of losses: " + str(cutoff[0]))
        print(
            "J == "
            + jump
            + " and L == "
            + str(cutoff[0])
            + ". Cutoff values are J => "
            + str(args.jumps)
            + " and L => "
            + str(args.losses)
        )
        print("No events detected, empty output file")