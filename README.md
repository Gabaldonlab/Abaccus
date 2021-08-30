# Abaccus v1.2

**Description**

Abaccus is a python3-based script for the detection of Horizontal Gene Transfer (HGT) events.

**Requirements:**

All dependencies should be installed before running the script: python (www.python.org) and ETE3 package (http://etetoolkit.org).


**Installation**

`git clone https://github.com/Gabaldonlab/Abaccus`

**Input:**

The program requires these main input:

* **-i**: a **phylogenetic tree** representing the evolution of a gene family. The **phylogenetic tree** can be provided in standard **newick** format, or as provided by the **phylomizer pipeline** (https://github.com/Gabaldonlab/phylomizer). In that case the input should be the folder of a phylomizer output. It is expected to include a phylogenetic tree in newick format used as the input, and a rank file with the score for the different models. It will select the best available newick tree based on the score.

* **-t**: The true distribution of species either as: a **.csv** file containing a **taxonomy** (i.e data/taxonomy.csv, included in this folder). The **taxonomy** should include at least all species contained in the **phylogenetic tree**. Or a species tree where the true phyloegenetic relationships between all the species present in the gene family tree are present. This file must end with **.newick**,**.nw** or **.nwk**.

* **-m**:
This is the name of the central sequence used for building the tree. By default, the script will assume that the name of the sequence is included in the name of the file containing the tree and will search it using a regular expression match. In any case, unless specified by the user, the program will assume that the name of the central sequence of the tree is contained in the name of the directory containing the tree (if used as input), and/or the newick tree analyzed. Currently, the script assumes that sequences use a Uniprot-like nomenclature system, with an alphanumerical ID value, followed by a mnemonic code of the species separated by underscore (I. e I7M603_TETTS; where I7MG03 is the Uniprot ID and TETTS is a mnemonic code for *Tetrahymena termophila*).

For the other arguments run `python scripts/abaccus.py -h`

**Output:**

Abaccus reports a raw text summarizing the number of loses needed to explain the phylogenetic distribution assuming exclusively vertical inheritance, the taxonomic levels at which the loses are predicted, and a printed version of the subtree containing the predicted transferred group as well as the branches in which it is embedded (sister branch and the sister branch to the common ancestor of the event + sister branch).

**The algorithm:**

TODO

**Examples:**

The current test directory contains three examples. Each directory contains the newick file, as well as the fasta file used for building the tree.

In order to launch the script with the examples, just use as input any of the trees in newick format. For example:
> python ./scripts/abaccus.py -i ./test/Example1/O59828_ALR1_SCHPO.tree.phyml.nj.LG.nw -t data/taxonomy.csv -o test/Example1/out_abaccus.txt

In this example, the program will create a file called abacus_output.txt with the results if it doesn't exists already. If it does, it will just append the output to the existing file. To send the results to STDOUT, use the -v or --verbose flag.

The included examples are the following:

  1) Example1 contains a well described horizontal gene transfer event in *Schizosaccharomyces pombe*, containing two highly related alanine racemases.

  2) Example2 contains a putative event described in Naranjo-Ortiz, et al (2016). Is a complex event involving a primary transfer of an aspartate-glutamate-hydantoin racemase from bacteria to fungi, and an uncertain number of secondary transferences between different fungal clades.

  3) Example3 contains a possible transference of an aspartate-glutamate-hydantoin racemase in the foraminiferan *Reticulomyxa
  filosa*. While the topology of the tree suggest a very possible horizontal gene transfer, abaccus will reject the event based on the lack of representation of related organisms. Negative control.

**Taxonomy file:**

The file *data/taxonomy.csv* is a manually curated taxonomy but the script *scripts/uniprot_parse.py* can build a NCBI based taxonomy starting from the gene family tree or from a file with a mnemonic code per line.

It depends on the file *data/speclist_24_08_21.txt* (if you want to update it simply run `curl https://www.uniprot.org/docs/speclist.txt -o data/speclist_$(date +'%d_%m_%y').txt`)

To use this script run
> python scripts/uniprot_parse.py -i $input_file -u da
ta/speclist_24_08_21.txt

By default this will create the file "taxonomy_mnemo.csv" in the present directory. The first time may take a while as ete3 needs to download the db.

By default the scripts fill missing data with the closest known more general phylum. This can be avoided with --nofill option. It also ignores prokaryotic entries as *Abaccus assumes that there are only eukaryotes in the taxonomy file*. This behaviour can be avoided with --proka. Finally, the script assumes a default scheme of ranks, which are: Species name, Genus, Family, Order, Class, Phylum, Kingdom, Superkingdom. However, with --max_cov option the script can detect which is the most common set of phyla and use that.
