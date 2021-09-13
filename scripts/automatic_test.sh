rm test/Example*/Ex*.txt

echo "Test 1 - R = 6, J = 9, L = 5"

python scripts/abaccus.py -i ./test/Example1/O59828_ALR1_SCHPO.tree.phyml.nj.LG.nw -t ./data/taxonomy.csv -o ./test/Example1/Example1_out.txt

grep "J =" ./test/Example1/Example1_out.txt

echo -e "\n"
echo "Test 2 - R = 6,J = 4,L = 6"

python scripts/abaccus.py -i ./test/Example2/B6Q7W7_B6Q7W7_PENMQ.tree.phyml.ml.LG.nw -t ./data/taxonomy.csv -o ./test/Example2/Example2_out.txt --image

grep "J =" ./test/Example2/Example2_out.txt

echo -e "\n"
echo "Test 3 - empty output"

python scripts/abaccus.py -i ./test/Example3/X6NCM0_X6NCM0_RETFI.tree.phyml.ml.LG.nw -t ./data/taxonomy.csv -o ./test/Example3/Example3_out.txt --image

echo -e "\n"
echo "Test 4 - Taxonomy"

echo "Creating taxonomy file"

python scripts/uniprot_parse.py -i ./test/Example_phylo/Phy0002LFP_CANAL.nwk -o ./test/Example_phylo/taxon_euka.csv -u data/speclist_24_08_21.txt

python scripts/abaccus.py -i ./test/Example_phylo/Phy0002LFP_CANAL.nwk -t ./test/Example_phylo/taxon_euka.csv -o ./test/Example_phylo/Example_tax.txt --image

grep "J =" ./test/Example_phylo/Example_tax.txt

echo -e "\n"
echo "Test 4 - Phylogeny mode, empty"

python scripts/abaccus.py -i ./test/Example_phylo/Phy0002LFP_CANAL.nwk -t ./test/Example_phylo/rooted_sprax_510.nwk -o ./test/Example_phylo/Example_phylo.txt -p test/Example_phylo/proka_mnemo.txt --image
