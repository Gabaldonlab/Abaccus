python scripts/abaccus.py -i ./test/Example1/O59828_ALR1_SCHPO.tree.phyml.nj.LG.nw -t ./data/taxonomy.csv -o ./test/Example1/Example1_out.txt


python scripts/abaccus.py -i ./test/Example2/B6Q7W7_B6Q7W7_PENMQ.tree.phyml.ml.LG.nw -t ./data/taxonomy.csv -o ./test/Example2/Example2_out.txt --image


python scripts/abaccus.py -i ./test/Example3/X6NCM0_X6NCM0_RETFI.tree.phyml.ml.LG.nw -t ./data/taxonomy.csv -o ./test/Example3/Example3_out.txt --image

# grep if j = 9 and l = 5 and other known results

# check if phylogeny mode works

python scripts/abaccus.py -i ./test/Example_phylo/Phy0002LFP_CANAL.nwk -t ./test/Example_phylo/taxon_euka.csv -o ./test/Example_phylo/Example_tax.txt --image


python scripts/abaccus.py -i ./test/Example_phylo/Phy0002LFP_CANAL.nwk -t ./test/Example_phylo/rooted_sprax_510.nwk -o ./test/Example_phylo/Example_phylo.txt -p test/Example_phylo/proka_mnemo.txt --image
