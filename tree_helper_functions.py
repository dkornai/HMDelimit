'''
THIS MODULE CONTAINS HELPER FUNCTIONS FOR DEALING WITH 
NEWICK AND ETE3 TREE DATA STRUCTURES.
'''
## DEPENDENCDIES
# STANDARD LIBRARY DEPENDENCIES
import copy
import warnings

# EXTERNAL LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

## TYPE HINTING 
from custom_types import Species_name
from custom_types import Tree_newick
from custom_types import Population_list
from custom_types import file_path

## SPECIALIZED HELPER FUNCTIONS
 
# takes a tree where only leaf nodes are named, and names all internal nodes

    # The naming of internal nodes is accomplished by combining the names of 
    # the descendant nodes. This is done from leaf to root. 

def name_Internal_nodes (
        tree:                   Tree
                        ) ->    Tree:

    for node in tree.traverse("postorder"):
        if len(node.name) == 0:
            newname = ""
            for count, descendants in enumerate(node.iter_descendants("levelorder")):
                if count == 0 or count == 1:
                    newname += descendants.name
            node.name = newname
    
    return tree

# returns a list of node pairs that are mergeable
def get_Mergeable_from_tree (
        tree:                       Tree, 
        pops:                       Population_list
                            ) ->    list[list[Species_name]]:

    tree = copy.deepcopy(tree)
    tree.prune(pops)

    mergepairs = []
    for node in tree.traverse("levelorder"):
        if len(list(node.iter_descendants("levelorder"))) == 2:
            pair = [leafnode.name for leafnode in node.iter_descendants("levelorder")]
            mergepairs.append(pair)
    
    return mergepairs

# returns the names of the node pairs that can be split
def get_Splitable_from_tree (
        tree:                       Tree, 
        pops:                       Population_list
                            ) ->    list[list[Species_name]]:

    splitpairs = []
    for node in tree.traverse("postorder"):
        if node.name in pops and len(list(node.iter_descendants("levelorder"))) != 0:
            pair = []
            for count, descendants in enumerate(node.iter_descendants("levelorder")):
                if count == 0 or count == 1:
                    pair.append(descendants.name)
            if set(pair).isdisjoint(pops):
                splitpairs.append(pair)
    
    return splitpairs

# small wrapper function that returns an ete3 tree in a newick formatted string
def tree_To_Newick  (
        tree:               Tree
                    ) ->    Tree_newick:

    return tree.write(format=9)

# write a newick tree to a file
def write_Tree  (
    tree:               Tree_newick,
    filename:           file_path,
                ):
    
    f = open(filename, "x")
    f.writelines([f"{tree}"])

# display the ASCII representation of the tree
def tree_ASCII  (
        tree:           Tree_newick
                ) ->    str:

    etetree = name_Internal_nodes(Tree(tree))

    return etetree.get_ascii()