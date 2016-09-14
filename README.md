# Abaccus
Abaccus pipeline for detection of HGT events. Version 1.0

Initial version of the script. The program assumes as input a folder containing a phylogenetic tree constructed by the pipeline phylomizer (several of them are included in the examples folder) and the location of a csv file containing the taxonomic information for the species (I.e taxonomy.csv, included in this folder). Currently, the script assumes that sequences use a Uniprot-like nomenclature system, with an alphanumerical ID value, followed by a mnemonic code of the species separated by underscore (I. e I7M603_TETTS; where I7MG03 is the Uniprot ID and TETTS is a mnemonic code for Tetrahymena termophila).
