# Abaccus
Abaccus is a python-based script for the detection of Horizontal Gene Transfer (HGT) events. Version 1.0

Initial version of the script. The publication of this script is meant to serve other researchers to reproduce the results of Naranjo-Ortiz et al (2016). We do not responsabilize of the problems that may arise if the script is used for other purposes.

Requires the Ete2 package (http://www.http://etetoolkit.org/)

The program assumes as input a folder containing a phylogenetic tree constructed by the pipeline phylomizer (https://github.com/Gabaldonlab/phylomizer) or a newick tree; and the location of a csv file containing the taxonomic information for the species (I.e taxonomy.csv, included in this folder, set as default). In any case, unless specified by the user, the program will assume that the name of the central sequence of the tree is contained in the name of the directory containing the tree (if used as input), and/or the newick tree analyzed. Currently, the script assumes that sequences use a Uniprot-like nomenclature system, with an alphanumerical ID value, followed by a mnemonic code of the species separated by underscore (I. e I7M603_TETTS; where I7MG03 is the Uniprot ID and TETTS is a mnemonic code for *Tetrahymena termophila*).

The output is a piece of raw text summarizing the number of loses needed to explain the phylogenetic distribution assuming exclusively vertical inheritance, the taxonomic levels at which the loses are predicted, and a printed version of the subtree containing the predicted transfered group as well as the branches in which it is embedded (sister branch and the sister branch to the common ancestor of the event + sister branch).

To see the rest of the options, just type:
> python ./abaccus.py -h





**Examples**

The current directory contains three examples. Each directory contains the newick file, as well as the fasta file used for building the tree.

In order to launch the script with the examples, just use as input any of the trees in newick format. For example:
> python ./abaccus.py -i ./Example1/O59828_ALR1_SCHPO.tree.phyml.nj.LG.nw

In this example, the program will create a file called abacus_output.txt with the results if it doesn't exists already. If it does, it will just append the output to the existing file. To send the results to STDOUT, use the -v or --verbose flag.

The included examples are the following:
  
  1) Example1 contains a well described horizontal gene transfer event in *Schizosaccharomyces pombe*, containing two highly related alanine racemases.
  
  2) Example2 contains a putative event described in Naranjo-Ortiz, et al (2016). Is a complex event involving a primary transfer of an aspartate-glutamate-hydantoin racemase from bacteria to fungi, and an uncertain number of secondary transferences between different fungal clades.
  
  3) Example3 contains a possible transference of an aspartate-glutamate-hydantoin racemase in the foraminiferan *Reticulomyxa 
  filosa*. While the topology of the tree suggest a very possible horizontal gene transfer, abaccus will reject the event based on the lack of representation of related organisms. Negative control.
