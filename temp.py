from samba_sampler import newick, tree2matrix, DistanceMatrix


# Test case
# TEST = "(((culi1244:0.216844,deni1241:0.216844):0.0787837,jama1261:0.29563):1.07742,paum1247:1.37305)"
# trees = newick.loads(TEST)
# tree = trees[0]
# distances = tree2matrix(tree)
# exit()


with open("src/samba_sampler/etc/global_tree.gled.newick", encoding="utf-8") as f:
    tree = newick.load(f)[0]

matrix = tree2matrix(tree)

matrix.save_matrix("gled.bz2")

# import samba_sampler as samba
# sampler = samba.MyLanguageSampler()
