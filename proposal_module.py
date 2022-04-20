'''
THIS MODULE CONTAINS SPECIALIZED FUNCTIONS FOR GENERATING NEW 
SPECIES DELIMITATION PROPOSALS. THE PROPOSALS CAN EITHER BE 
TO SPLIT POPULATIONS, OR TO MERGE POPULATIONS.
'''

## DEPENDENCDIES
# STANDARD LIBRARY DEPENDENCIES
import copy
import warnings

# EXTERNAL LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

# HELPER FUNCTION DEPENDENCIES
from helper_functions import flatten

## TYPE HINTING 
from custom_types import Alias_dict
from custom_types import Species_name
from custom_types import Tree_newick
from custom_types import Population_list
from custom_types import HM_mode
from custom_types import Imap_list


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

# small wrapper function that returns an ete3 tree in a newick formatted string
def tree_To_Newick  (
        tree:               Tree
                    ) ->    Tree_newick:

    return tree.write(format=9)

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

# remaps the leaf nodes of the tree to the most recent ancestor that is currently in the pop list
def remap_from_tree (
        tree:               Tree, 
        pops:               Population_list
                    ) ->    Alias_dict:

    remap_dict = {}
    for node in tree.traverse("levelorder"):
        # for each leaf
        if node.is_leaf() == True:
            p_labels = []
            p_labels.append(node.name)
            for ancestor in node.get_ancestors():
                p_labels.append(ancestor.name)
            # find the first ancestor with a name in the pops list
            for x in p_labels:
                for y in list(reversed(pops)):
                    if x == y:
                        remap_dict[p_labels[0]] = y
                        break 
                else:
                    continue
                break
    
    return remap_dict

# writes the remap proposed by the program to a new imap list
def remap_to_imapList   (
        imap_base:              Imap_list, 
        remap:                  Alias_dict
                        ) ->    Imap_list:

    IDs = []
    pop = []
    for ind in imap_base:
        IDs.append(ind)
        pop.append(remap[imap_base[ind]])

    return [IDs, pop]


## MAIN PROPOSAL FUNCTIONS

# generate the starting state of accepted populations for the HM stage, depending on if the method is progressive
# merge, or progressive split. Also return the number of populations when the end condition is reached
def get_HM_StartingState(
        input_guide_tree:       Tree_newick, 
        mode:                   HM_mode
                        ) ->    tuple[Population_list, int]:
    
    tree = Tree(input_guide_tree)
    tree = name_Internal_nodes(tree)
    

    if mode == "merge":
        # the starting configuration is that all nodes are accepted as species, and then progressively rejected
        starting_pops = [node.name for node in tree.traverse("preorder")]
        halt_pop_number = 1 # the program should always end if all all nodes except the root have been rejeceted
    elif mode == "split":
        # the starting configurateion is that all nodes except for the root are rejected as species
        starting_pops = [node.name for x, node in enumerate(tree.traverse("preorder")) if x == 0]
        halt_pop_number = len(tree) # the program should always end if all possible nodes were split
    
    return starting_pops, halt_pop_number

# print(get_HM_StartingState("(((A, B), CD), X);", "merge"))
# print(get_HM_StartingState("(((A, B), CD), X);", "split"))


# generate the proposal for which populations should be removed or added to the accepted list, depending on the mode

    # This function is used in the Hierarchical Method to generate the new and topology and IMAP,
    # corresponding to the next iteration of the process.
    
def HMproposal  (
        guide_tree_newick:  Tree_newick, 
        base_indpop_dict, 
        current_pops_list:  Population_list, 
        mode:               HM_mode
                ) ->        tuple[list[list[Species_name]], Tree_newick, Imap_list]:

    # transform newick tree to ete3, and name all the internal nodes
    guide_ete3_tree = Tree(guide_tree_newick)
    guide_ete3_tree = name_Internal_nodes(guide_ete3_tree)

    # In merge mode, the Tau and Theta values are calculated using the current topology.
    if mode == "merge":
        prop_change = get_Mergeable_from_tree(guide_ete3_tree, current_pops_list)
        proposed_pops = current_pops_list

    # In split mode, Tau and Theta values are calculated using a proposed split topology.    
    elif mode == "split":
        prop_change = get_Splitable_from_tree(guide_ete3_tree, current_pops_list)
        proposed_pops = current_pops_list + flatten(prop_change)
    
    proposed_remap = remap_from_tree(guide_ete3_tree, proposed_pops)
    
    # imap corresponding to the new proposal
    proposed_imap = remap_to_imapList(base_indpop_dict, proposed_remap)
    
    # tree topology corresponding to the new proposal
        # This is constucted by pruning the guide tree to only include the proposed populations
    proposed_tree = copy.deepcopy(guide_ete3_tree)
    new_tree = proposed_tree
    new_tree.prune(proposed_pops)
    proposed_tree = tree_To_Newick(new_tree)

    # return the three components of the proposal
    return prop_change, proposed_tree, proposed_imap