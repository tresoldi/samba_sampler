#!/usr/bin/env python3

"""
Build supplementary trees.
"""

# Import Python libraries
from ete3 import Tree
import glob
import csv
from pathlib import Path

BASE_PATH = Path(__file__).parent.parent

# list of all tree files
tree_files = glob.glob('gled.trees/*.tree')

max_branch_length = 0

# iterate over all tree files to find the maximum branch length
for tree_file in tree_files:
    t = Tree(str(tree_file))
    for node in t.traverse("postorder"):
        if node.get_distance(t) > max_branch_length:
            max_branch_length = node.get_distance(t)

print("max branch length:", max_branch_length)

# Read the glottolog.tsv file
glottolog_dump = sorted(BASE_PATH.glob('src/samba_sampler/etc/glottolog.*.tsv'))[-1]
isolates = []
with open(glottolog_dump, 'r', encoding="utf-8") as f:
    tsv_reader = csv.reader(f, delimiter='\t')
    headers = next(tsv_reader)  # get the headers
    for row in tsv_reader:
        row_dict = dict(zip(headers, row))
        if row_dict['isolate'] == 'True':
            isolates.append(row_dict['glottocode'])

# create the global tree
global_tree = Tree()

# add each tree from disk to the global tree
for tree_file in tree_files:
    t = Tree(tree_file)
    # Change the leaf labels to just the Glottocode
    for leaf in t.iter_leaves():
        leaf.name = leaf.name.split("_")[1]
    t.dist = max_branch_length
    global_tree.add_child(t)

# add isolates to the global tree
for isolate in isolates:
    t = Tree(name=isolate)
    t.dist = max_branch_length
    global_tree.add_child(t)

# print the global tree
#print(global_tree)

# count the number of leaves
num_taxa = len(global_tree.get_leaves())

print(f'The number of taxa in the tree is: {num_taxa}')

# write the global tree to a file
output = BASE_PATH / "src" / "samba_sampler" / "etc" / "global_tree.gled.newick"
global_tree.write(outfile=output)