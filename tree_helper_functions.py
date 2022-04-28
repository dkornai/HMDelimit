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
    from ete3 import TreeStyle
    from ete3 import NodeStyle
    from ete3 import TextFace
    from ete3 import AttrFace

# HELPER DEPENDENCIES
from helper_functions import flatten
from helper_functions import extract_Name_TauTheta_dict

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

## PROVIDE A VISUALIZATION OF THE CHANGES THAT WERE MADE IN THE HM STAGE
'''
This section is used to provide a visual output of the decisions made during the HM stage.
The requested tree is visualized, and the populations that would be changed are colored
according to whether they were accepted to change or not.
'''
# helper datasets for styling the feedback tree
general_style = NodeStyle()
general_style["size"] = 0
general_style["hz_line_width"] = 2
general_style["vt_line_width"] = 2

accepted_leaf = NodeStyle()
accepted_leaf["size"] = 0
accepted_leaf["hz_line_color"] = "LimeGreen"
accepted_leaf["hz_line_width"] = 4
accepted_leaf["hz_line_type"] = 1
accepted_leaf["vt_line_width"] = 2

accepted_node = NodeStyle()
accepted_node["size"] = 0
accepted_node["hz_line_width"] = 2
accepted_node["vt_line_color"] = "LimeGreen"
accepted_node["vt_line_width"] = 4
accepted_node["vt_line_type"] = 1

rejected_leaf = NodeStyle()
rejected_leaf["size"] = 0
rejected_leaf["hz_line_color"] = "Grey"
rejected_leaf["hz_line_width"] = 4
rejected_leaf["hz_line_type"] = 0
rejected_leaf["vt_line_width"] = 2

rejected_node = NodeStyle()
rejected_node["size"] = 0
rejected_node["vt_line_color"] = "Grey"
rejected_node["vt_line_width"] = 4
rejected_node["vt_line_type"] = 0
rejected_node["hz_line_width"] = 2

ts = TreeStyle()
ts.branch_vertical_margin = 30
ts.show_scale = False
ts.margin_bottom = 10
ts.margin_left = 10
ts.margin_right = 10
ts.margin_top = 10
ts.scale = 100

# main function implementing the drawing of the feedback tree
def visualize_decision(proposed_tree, MSC_param, BPP_outfile, proposed_changes, decision):
    
    # collect the accepted and rejected changes
    if len(decision) == 0:
        accepted_changes = []
    else:
        accepted_changes = flatten(decision)
    all_changes = flatten(proposed_changes)
    rejected_changes = []
    for item in all_changes:
        if item not in accepted_changes:
            rejected_changes.append(item)

    # collect the split ages
    tau_dict = extract_Name_TauTheta_dict(BPP_outfile)[0]
    theta_dict = extract_Name_TauTheta_dict(BPP_outfile)[1]

    # collect the decision parameters that are going to be visualized
    gdi_dict = {}
    age_dict = {}
    for pair in MSC_param:
        pairname = str(pair)[2:-2].split("', '")
        merged_name = "".join(str(pair)[2:-2].split("', '"))
        if MSC_param[pair]["gdi_1"] != "?":
            gdi_dict[pairname[0]] = MSC_param[pair]["gdi_1"]
        if MSC_param[pair]["gdi_2"] != "?":
            gdi_dict[pairname[1]] = MSC_param[pair]["gdi_2"]
        if MSC_param[pair]["age"] != "?":
            age_dict[merged_name] = MSC_param[pair]["age"]

    # prepare the tree for visualization
    tree = Tree(proposed_tree)
    tree = name_Internal_nodes(tree)

    # put the data onto nodes
    for i, node in enumerate(tree.traverse("preorder")):
        name = node.name
        if name in gdi_dict:
            node.add_features(gdi=f" GDI={gdi_dict[name]}")
        if name in tau_dict:
            node.add_features(tauval=f" τ={tau_dict[name]}")
        if name in theta_dict:
            node.add_features(theta=f"θ={theta_dict[name]}")
        if name in age_dict:
            node.add_features(age=f" age={age_dict[name]}")

    # make the tree ultrametric
    for node in tree.traverse("levelorder"):
        if node.name in tau_dict:
            node.add_feature(pr_name="tau", pr_value=tau_dict[node.name])

    tau_sums = []
    for node in tree.traverse():
        if node.is_leaf():
            ancestors = node.get_ancestors()
            tau_sum = 0
            for ancestor in ancestors:
                tau_sum += ancestor.tau
            tau_sums.append(tau_sum)

    multiplier = (1/max(tau_sums))

    for node in tree.iter_descendants("postorder"):
        ancestor = node.up
        ancestor_name = ancestor.name
        node.dist = (tau_dict[ancestor_name]*multiplier*10)

    target_dist = []
    for node in tree.traverse():
        if node.is_leaf():
            ancestors = node.get_ancestors()
            dist_sum = node.dist
            for ancestor in ancestors:
                dist_sum += ancestor.dist
            target_dist.append(dist_sum)
    
    target = max(target_dist)

    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            node_dist_to_root = 0
            ancestors = node.get_ancestors()
            for ancestor in ancestors:
                node_dist_to_root += ancestor.dist
            node.dist = (target-node_dist_to_root)

    # scale theta values to be well visible
    leaf_names = [node.name for node in tree.traverse() if node.is_leaf()]
    max_leaf_theta = 0
    for item in leaf_names:
        if theta_dict[item] > max_leaf_theta:
            max_leaf_theta = theta_dict[item]

    theta_mult = 50/max_leaf_theta # scale so that the largest theta corresponds to a radius of 50

    theta_dict = {key:int(theta_dict[key]*theta_mult) for key in theta_dict}
    for key in theta_dict:
        if theta_dict[key] < 1:
            theta_dict[key] = 1

    # get directory of left and right nodes, this is needed to place annotations correctly
    left_nodes = []
    right_nodes = []
    for i, node in enumerate(tree.traverse("levelorder")):
        if i != 0 and i%2 == 1:
            left_nodes.append(node.name)
        if i != 0 and i%2 == 0:
            right_nodes.append(node.name)

    # style the tree with the base style
    for node in tree.traverse():
        node.set_style(general_style)

    # style the tree to reflect the decision
    age_ancestors = []
    tau_ancestors = []
    
    for node in tree.traverse("preorder"):
        name = node.name
        if node.is_leaf():
            
            # write in GDI value if available
            if name in gdi_dict:
                face = AttrFace("gdi", fsize=10, fstyle="italic")
                face.hz_align = 0
                node.add_face(face, column=1, position = "aligned")

            ancestor = node.up
            ancestorname = ancestor.name
            # if a node is the ancestor to a leaf node pair
            if len(list(ancestor.iter_descendants("levelorder"))) == 2:
                # write in age in generations
                if ancestorname in age_dict and ancestorname not in age_ancestors:
                    ageface = AttrFace("age", fsize=10, fstyle="italic")
                    ageface.margin_right = -10*(len(str(age_dict[ancestorname]))+5)
                    ageface.vt_align = 1
                    ancestor.add_face(ageface, column=2, position = "branch-bottom")
                    age_ancestors.append(ancestorname)
                # write in tau value
                if ancestorname in tau_dict and ancestorname not in tau_ancestors:
                    tauface = AttrFace("tauval", fsize=10)
                    tauface.margin_right = -90
                    tauface.vt_align = 1
                    ancestor.add_face(tauface, column=2, position = "branch-top")
                    tau_ancestors.append(ancestorname)

            # edit branches and circle to accepted style
            if name in accepted_changes:
                node.set_style(accepted_leaf)
                ancestor = node.up
                ancestor.set_style(accepted_node)

                if name in theta_dict:
                    nstyle = NodeStyle()
                    nstyle["hz_line_width"] = 4
                    nstyle["hz_line_color"] = "LimeGreen"
                    nstyle["vt_line_color"] = "LimeGreen"
                    nstyle["vt_line_width"] = 4
                    nstyle["vt_line_type"] = 1
                    nstyle["hz_line_type"] = 1
                    nstyle["fgcolor"] = "LimeGreen"
                    nstyle["shape"] = "sphere"
                    nstyle["size"] = theta_dict[name]
                    node.set_style(nstyle)

            # edit branches and circle to rejected style
            elif name in rejected_changes:
                node.set_style(rejected_leaf)
                ancestor = node.up
                ancestor.set_style(rejected_node)
                if name in theta_dict:
                    nstyle = NodeStyle()
                    nstyle["hz_line_color"] = "Grey"
                    nstyle["hz_line_width"] = 4
                    nstyle["hz_line_type"] = 0
                    nstyle["vt_line_width"] = 2
                    nstyle["fgcolor"] = "Grey"
                    nstyle["shape"] = "sphere"
                    nstyle["size"] = theta_dict[name]
                    node.set_style(nstyle)
            
            # add in theta circle for remaining reaves
            elif name in theta_dict:
                nstyle = NodeStyle()
                nstyle["hz_line_width"] = 2
                nstyle["vt_line_width"] = 2
                nstyle["fgcolor"] = "Gainsboro"
                nstyle["shape"] = "circle"
                nstyle["size"] = theta_dict[name]
                node.set_style(nstyle)
            
            # add in actual theta value
            if name in theta_dict:
                thetaface = AttrFace("theta", fsize=10)
                node.add_face(thetaface, column=1, position = "branch-top")
            
    tree.render("Decision_Visualization.png", tree_style=ts)