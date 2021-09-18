# Abaccus v1.2

**Description**

Abaccus is a python3-based script for the detection of Horizontal Gene Transfer (HGT) events.

**Requirements:**

All dependencies should be installed before running the script: python (www.python.org) and ETE3 package (http://etetoolkit.org).


**Installation**

`git clone https://github.com/Gabaldonlab/Abaccus`

**Input:**

The program requires these main input:

* **-i**: a **phylogenetic tree** representing the evolution of a gene family. The **phylogenetic tree** can be provided in standard **newick** format, or as provided by the [**phylomizer pipeline**](https://github.com/Gabaldonlab/phylomizer). In that case the input should be the folder of a phylomizer output. It is expected to include a phylogenetic tree in newick format used as the input, and a rank file with the score for the different models. It will select the best available newick tree based on the score.

* **-n**: the naming scheme of the gene names. Available options are pdb (phylomeDB, Default, gene_species) or uniref (species after the last underscore).

* **-t**: The true distribution of species either as: a **.csv** file containing a **taxonomy** (i.e data/taxonomy.csv, included in this folder). The **taxonomy** should include at least all species contained in the **phylogenetic tree**. Or a species tree where the true phyloegenetic relationships between all the species present in the gene family tree are present. This file must end with **.newick**, **.nw** or **.nwk**.

* **-m**:
This is the name of the central sequence used for building the tree. By default, the script will assume that the name of the sequence is included in the name of the file containing the tree and will search it using a regular expression match. In any case, unless specified by the user, the program will assume that the name of the central sequence of the tree is contained in the name of the directory containing the tree (if used as input), and/or the newick tree analyzed. Currently, the script assumes that sequences use a Uniprot-like nomenclature system, with an alphanumerical ID value, followed by a mnemonic code of the species separated by underscore (I. e I7M603_TETTS; where I7MG03 is the Uniprot ID and TETTS is a mnemonic code for *Tetrahymena termophila*).

* **-r**:
This is the number of round Abaccus will do before exiting. Meaning how much deep inside the tree it will search for HGT events. The Default is 6. If this number is set too high and Abaccus eventually reaches the root the actual rounds will be less.

* **--collapse_dup**: Collapse lineage specific expansions in the genetree (it will be done randomly but the main sequence is guaranteed to be included).

For the other arguments run `python scripts/abaccus.py -h`

**Output:**

Abaccus reports a raw text summarizing the number of loses needed to explain the phylogenetic distribution assuming exclusively vertical inheritance, the taxonomic levels at which the loses are predicted, and a printed version of the subtree containing the predicted transferred group as well as the branches in which it is embedded (sister branch and the sister branch to the common ancestor of the event + sister branch).

If `--image` option is used and Abaccus finds an event an image of the tree will be produced in output directory called outfile_abaccus.png. The seed will have a green dot, the blue node is the start point and the red dot is the end point.

**The algorithm:**

###### Taxonomy mode:

"We first root the tree at the farthest leaf from the seed
protein found in the tree. Then we run through every node
from the seed protein to the root node (depending on Round parameter). For each node we
determine the taxonomic classification of the node and its parent
node. Then we compare the two taxonomic classifications. The
“jump” parameter (J) is defined as the difference between the
taxonomic level found at the parent node and the one found
at the current node. As seen in Figure 1A the jump parameter
between F. oxysporum and F. graminearum is equal to 1 because
we move from species level to genus level while in Figure 1B
the jump parameter between Fusarium and A. nidulans is equal
to four because we jump from genus level to phylum level. We
then compute the minimal number of loss events (L) between the
node and its sister node. We use a very parsimonious approach
that infers that for each taxonomic level of difference between a
node and its parental node one single loss event has happened
only if there is at least a species in our database belonging to
that taxonomic classification level that is not present in the nodes.
This also implies that no loss is inferred if no other member of a
given taxonomic category is present in the database."

###### Phylogeny mode:

Initially, we compute the "species to age" dictionary starting from the seed species.
This dictionary contains information about the distance between the seed
species and the MRCA + 1 with all other species.
The gene tree is then rooted with this dictionary. This means that the root will be the furthest and oldest gene from the furthest species.

Then, starting from the species tree, we consider internal nodes as proxies of taxonomic groupings. For each node we compare the phylogenetic classification of the node and its parent node. The “jump” parameter (J) is defined as maximum value of the species2age dictionary between all species in the analyzed node.
If the threshold for J is reached Losses will be computed this way:
Traverse the previously found suspect gene sub-tree and get all nodes comprising the seed sequence. For each of this nodes a loss will be added if the species in the gene tree node will not be all the species present in the corresponding species tree node.


In both cases, the algorithm will run through every node either until it finds an event or until it makes N rounds (N=6 by default). If N is set to 6 and the algorithm reaches the tree root before these round are over the algorithm will stop.

**Examples:**

The current test directory contains three examples. Each directory contains the newick file, as well as the fasta file used for building the tree.

In order to launch the script with the examples, just use as input any of the trees in newick format. For example:
> python ./scripts/abaccus.py -i ./test/Example1/O59828_ALR1_SCHPO.tree.phyml.nj.LG.nw -t data/taxonomy.csv -o test/Example1/out_abaccus.txt

In this example, the program will create a file called abacus_output.txt with the results if it doesn't exists already. If it does, it will just append the output to the existing file. To send the results to STDOUT, use the -v or --verbose flag.

or simply do to run all examples:

`./scripts/automatic_test.sh` to launch all 4 examples.


The included examples are the following:

  1) Example1 contains a well described horizontal gene transfer event in *Schizosaccharomyces pombe*, containing two highly related alanine racemases.

  2) Example2 contains a putative event described in Naranjo-Ortiz, et al (2016). Is a complex event involving a primary transfer of an aspartate-glutamate-hydantoin racemase from bacteria to fungi, and an uncertain number of secondary transferences between different fungal clades.

  3) Example3 contains a possible transference of an aspartate-glutamate-hydantoin racemase in the foraminiferan *Reticulomyxa
  filosa*. While the topology of the tree suggest a very possible horizontal gene transfer, abaccus will reject the event based on the lack of representation of related organisms. Negative control.

  4) Example with a gene tree from Phylome 501. Both with the automatic taxonomy and with the species tree.

**Taxonomy file:**

The file *data/taxonomy.csv* is a manually curated taxonomy but the script *scripts/uniprot_parse.py* can build a NCBI based taxonomy starting from the gene family tree or from a file with a mnemonic code per line.

It depends on the file *data/speclist_24_08_21.txt* (if you want to update it simply run `curl https://www.uniprot.org/docs/speclist.txt -o data/speclist_$(date +'%d_%m_%y').txt`)

To use this script run
> python scripts/uniprot_parse.py -i $input_file -u data/speclist_24_08_21.txt

By default this will create the file "taxonomy_mnemo.csv" in the present directory. The first time may take a while as ete3 needs to download the db.

By default the scripts fill missing data with the closest known more general phylum. This can be avoided with --nofill option. It also ignores prokaryotic entries as *Abaccus assumes that there are only eukaryotes in the taxonomy file*. This behaviour can be avoided with --proka. Finally, the script assumes a default scheme of ranks, which are: Species name, Genus, Family, Order, Class, Phylum, Kingdom, Superkingdom. However, with --max_cov option the script can detect which is the most common set of phyla and use that.
