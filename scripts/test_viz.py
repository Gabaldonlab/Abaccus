from ete3 import random_color, Tree, faces, TreeStyle
import sys
sys.path.insert(1, "scripts")
import functions as fn

taxofile = "test/phylo_data/proka_tax.csv"
taxo_dict = fn.taxonomist(taxofile)
new_dict = {}
for key in taxo_dict:
    species = taxo_dict[key][0]
    rest = taxo_dict[key][1:]
    new_dict[species] = rest

ref_tree = "test/phylo_data/rooted_sprax_510.nwk"
sptree = Tree(ref_tree)

add_set = set()
for el in new_dict:
    for a in new_dict[el]:
        add_set.add(a)

color_dict = {}
colors = random_color(num=len(set(add_set)))
color_dict = {el[1]:el[0] for el in list(zip(colors, list(add_set)))}

for key in new_dict:
    new_tuple = [(el, color_dict[el]) for el in new_dict[key]]
    new_dict[key] = new_tuple

for node in sptree.iter_leaves():
    if node.name in taxo_dict.keys():
        node.species = taxo_dict[node.name][0]
        sp_df = new_dict[node.species]
        node.add_feature("col_df", sp_df)


def layout_species(node):
    width = 100 # try to find a mulitplicator or something
        # If node is a leaf, add the nodes name and a its scientific
        # name
    node.img_style["size"] = 0
    if node.is_leaf():
        name_face = faces.AttrFace("species")
        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
        col_idx = 1
        for clade in node.col_df:
            rect = faces.RectFace(width, 20, bgcolor=clade[1], fgcolor=clade[1], label={"text":clade[0], "color":"white", "fontsize":6})
            faces.add_face_to_node(rect, node, column=col_idx, aligned=True)
            col_idx += 1

ts = TreeStyle()
ts.show_leaf_name = False
ts.allow_face_overlap = True
ts.draw_aligned_faces_as_table = True
ts.layout_fn = layout_species
#sptree.render("prova.png",tree_style = ts)
sptree.show(tree_style = ts)


input = "test/phylo_data/Phy0002LFP_CANAL.nwk"
central_seq = "Phy0002LFP_CANAL"

phylotree = Tree(input)

def layout_genes(node):
        # If node is a leaf, add the nodes name and a its scientific
        # name
    if node.is_leaf():
        nameFace = faces.TextFace(node.name)
        faces.add_face_to_node(nameFace, node, column=0)
        if node.name == central_seq:
            node.img_style["size"] = 7
            node.img_style["shape"] = "circle"
            node.img_style["fgcolor"] = "darkred"
        else:
            node.img_style["size"] = 0

ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = True
ts.show_branch_support = True
ts.layout_fn = layout_genes
phylotree.show(tree_style = ts)
