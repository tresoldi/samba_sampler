from samba_sampler import newick, tree2matrix, DistanceMatrix


# Test case
TEST = "(((culi1244:0.216844,deni1241:0.216844):0.0787837,jama1261:0.29563):1.07742,paum1247:1.37305)"
trees = newick.loads(TEST)
tree = trees[0]

distances = tree2matrix(tree)

exit()



with open("src/samba_sampler/etc/global_tree.gled.newick", encoding="utf-8") as f:
    tree = newick.load(f)[0]

lects = [node.name for node in tree.get_leaves()]
distances = tree2matrix(tree)

m = DistanceMatrix(lects)
for key, value in distances.items():
    print("***", key, value)
    lang1, lang2 = key
    print ("*", lang1, lang2, value)
    m.set(lang1, lang2, value)

m.save_matrix("gled.bz2")

print( m.keys[:20]  )

exit()



print(tree)

distances = newick2matrix(tree)
print("dist", distances)

exit()


print(culina.ancestors)
print(jama.ancestors)


#import samba_sampler as samba

#sampler = samba.MyLanguageSampler()