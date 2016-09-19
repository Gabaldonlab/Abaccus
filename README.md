# Abaccus
Abaccus pipeline for detection of HGT events. Version 1.0

Initial version of the script. 

The program assumes as input a folder containing a phylogenetic tree constructed by the pipeline phylomizer (several of them are included in the examples folder) and the location of a csv file containing the taxonomic information for the species (I.e taxonomy.csv, included in this folder). Currently, the script assumes that sequences use a Uniprot-like nomenclature system, with an alphanumerical ID value, followed by a mnemonic code of the species separated by underscore (I. e I7M603_TETTS; where I7MG03 is the Uniprot ID and TETTS is a mnemonic code for *Tetrahymena termophila*).

The output is a text summarizing the number of loses needed to explain the phylogenetic distribution assuming exclusively vertical inheritance, the taxonomic levels at which the loses are predicted, and a printed version of the subtree containing the predicted transfered group as well as the branches in which it is embedded (sister branch and the sister branch to the common ancestor of the event + sister branch).

In order to launch the script with the examples, just use as input any of the trees in newick format. For example:
> python ./abaccus.py -i ./Example1/O59828_ALR1_SCHPO.tree.phyml.nj.LG.nw

In this example, the program will create a file called abacus_output.txt with the results if it doesn't exists already. If it does, it will just append the output to the existing file. 

To see the rest of the options, just type:
> python ./abaccus.py -h
