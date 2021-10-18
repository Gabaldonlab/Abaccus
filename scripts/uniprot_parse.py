from ete3 import NCBITaxa, Tree
import pandas as pd
import argparse
import json

# TODOS:
# nicer implementation of autofill and is the naive one good?
# bacteria rarely have kingdom so maybe manual check regarding that
# add this to readme "curl https://www.uniprot.org/docs/speclist.txt -O"
# to download speclist
# works with ncbitaxa!

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default=False,
        required=True,
        help="The gene family tree or a txt file containing the mnemonic codes",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="./taxonomy_mnemo.csv",
        help="Name of the output file.",
    )
    parser.add_argument(
        "-u",
        "--uniprot_df",
        # required=True,
        help="File of the uniprot translation codes: https://www.uniprot.org/docs/speclist.txt",
    )
    parser.add_argument(
        "--nofill",
        default=False,
        action="store_true",
        help="Flag to avoid automatic filling",
    )
    parser.add_argument(
        "--proka",
        default=False,
        action="store_true",
        help="Flag to include also prokaryotic entries in output csv",
    )
    parser.add_argument(
        "--max_cov",
        default=False,
        action="store_true",
        help="Flag to use most represented set of taxonomic entries instead of default taxonomic entities",
    )

ncbi = NCBITaxa()
args = parser.parse_args()

input = args.input

taxo_dict = {}

already_done = False

if input.endswith(".txt"):
    with open(input) as inp:
        sp = [line.strip() for line in inp.readlines()]
        # remove duplicates
        sp = list(set(sp))
elif input.endswith(("nwk", "nw", "newick")):
    tree = Tree(input)
    mnemoset = {str(leaf)[str(leaf).rfind("_") + 1 :] for leaf in tree.get_leaf_names()}
    sp = [element for element in mnemoset]
elif input.endswith("json"):
    with open(input) as td:
        taxo_dict = json.load(td)
        taxo_dict = {v: k for k, v in taxo_dict.items()}
        already_done = True

if not already_done:
    sp_str = [el for el in sp if not el.isdecimal()]
    sp_num = [el for el in sp if el.isdecimal()]

    with open(args.uniprot_df) as f:
        targets = [line for line in f for s in sp_str if s in line]
        for line in targets:
            line = line.split()
            if len(line) > 2:
                mnemo = line[0]
                taxid = line[2].replace(":", "")
                if mnemo in sp:
                    taxo_dict[taxid] = mnemo
                elif taxid in sp:
                    taxo_dict[taxid] = mnemo

    for el in sp_num:
        if len(ncbi.get_taxid_translator([el])) > 0:
            taxo_dict[el] = el

    absent = set([abs for abs in sp if abs not in taxo_dict.values()])

    if absent:
        print("Warning: " + " ".join(absent) + " was not found in the dictionary.")
        print("You may want to add it manually")
    abs = [el for el in sp if el not in list(taxo_dict.values()) and not el.isdecimal()]
    abs_tx = [el for el in sp if el not in list(taxo_dict.keys()) and el.isdecimal()]

    abs = abs + abs_tx

    if len(abs) > 0:
        print(
            "Warning: "
            + str(len(abs))
            + " sequences were not found in "
            + str(args.uniprot_df)
        )
        print(abs)
        print(
            "Check if these are eukaryotes and manually annotate the csv, otherwise results may be misleading!"
        )

if not args.max_cov:
    # txid;mnemo;species;genre;family;Order;Class;Phylum;Kingdom;Superkingdom

    cols = [
        "#0-TaxaID",
        "#1-Mnemonic",
        "#2-Species name",
        "#3-Genus",
        "#4-Family",
        "#5-Order",
        "#6-Class",
        "#7-Phylum",
        "#8-Kingdom",
        "#9-Superkingdom",
    ]

    data = []

    for key, value in taxo_dict.items():
        txid = key
        mnemo = value
        try:
            sp_name = "".join(ncbi.get_taxid_translator([key]).values())
            lineage = ncbi.get_taxid_translator(ncbi.get_lineage(key))
            rank = ncbi.get_rank(lineage)

            final_dict = {}
            for k in lineage.keys():
                new_k = rank[k]
                final_dict[new_k] = lineage[k]
            supk = final_dict.get("superkingdom")

            # if supk == 'Eukaryota':
            sp = final_dict.get("species")
            gen = final_dict.get("genus")
            fam = final_dict.get("family")
            ord = final_dict.get("order")
            cla = final_dict.get("class")
            phy = final_dict.get("phylum")
            kin = final_dict.get("kingdom")

            row = [txid, mnemo, sp, gen, fam, ord, cla, phy, kin, supk]

            if args.proka:
                data.append(row)
            elif supk == "Eukaryota":
                data.append(row)

        except ValueError:
            print("Warning:" + str(txid) + "was not found")

    # warning there re n none
    df = pd.DataFrame(data=data, columns=cols)
else:
    # this is the script to do the longer taxonomy (VERY UGLY, to fix)
    data = []
    listone = []
    for key, value in taxo_dict.items():
        sp_name = "".join(ncbi.get_taxid_translator([key]).values())
        ordered = ncbi.get_lineage(key)
        lineage = ncbi.get_taxid_translator(ordered)
        rank = ncbi.get_rank(lineage)
        final_dict = []
        if lineage.get(2759):  # eukaryotes, hardcoded = bad
            for k in ordered:
                new_k = rank[k]
                final_dict.append((new_k, lineage[k]))
            final_dict.append(("mnemo", value))
            final_dict.append(("taxaid", key))
            listone.append(final_dict)
        else:
            next
    # delete clade and no rank which are not good for automization
    for el in listone:
        for val in list(el):
            if val[0] == "clade" or val[0] == "no rank":
                el.remove(val)
    # now we find the most common set of attributes
    # this may be simplified a lot!
    list_sets = []
    for el in listone:
        el_set = []
        for i in el:
            el_set.append(i[0])
        list_sets.append(el_set)
    used = []
    unique = [used.append(x) for x in list_sets if x not in used]
    max_occ = max([list_sets.count(el) for el in used])
    cols_set = [el for el in used if list_sets.count(el) == max_occ][0]

    for el in listone:
        el_dict = {i[0]: i[1] for i in el}
        row = []
        for i in cols_set:
            cell = el_dict.get(i)
            row.append(cell)
        data.append(row)

    df = pd.DataFrame(data=data, columns=cols_set)
    df = df[df.columns[::-1]]
    df = df.replace(r"\\n", " ", regex=True)


# loop over columns and find na
if args.nofill:
    print("There are these missing data!")
    print(df.isnull().sum())
else:
    print("There are these missing data!")
    print(df.isnull().sum())
    print("Filling missing data with the first nearest (less specific) clade")
    nas = df.isnull().sum().sum()
    while nas > 0:
        for col in range(2, df.shape[1]):
            row = df.index[df[df.columns[col]].isnull()].tolist()
            for el in row:
                # should it take the clade before or after
                # maybe to be the most conservative the one before
                # but I have to think about it
                # anyway a nicer thing could be implemented
                df.iloc[el, col] = df.iloc[el, col + 1]
        nas = df.isnull().sum().sum()

df.to_csv(args.output, index=False, header=True, sep=";")
