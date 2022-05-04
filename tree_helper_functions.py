'''
THIS MODULE CONTAINS HELPER FUNCTIONS FOR DEALING WITH 
NEWICK AND ETE3 TREE DATA STRUCTURES.
'''
## DEPENDENCDIES
# STANDARD LIBRARY DEPENDENCIES

import copy
import warnings
# set qtl to nonwindowed mode, this way the pipeline should work through the command line
import os
os.environ['QT_QPA_PLATFORM']='offscreen'

# EXTERNAL LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree
    from ete3 import TreeStyle
    from ete3 import NodeStyle
    from ete3 import TextFace
    from ete3 import AttrFace

import distinctipy

# HELPER DEPENDENCIES
from helper_functions import flatten
from helper_functions import string_limit
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
        tree:           Tree_newick,
        filename:       file_path,
                ):
    
    f = open(filename, "x")
    f.writelines([f"{tree}"])

# display the ASCII representation of the tree
def tree_ASCII  (
        tree:           Tree_newick
                ) ->    str:

    # start by naming the internal nodes, and assuming they are reasonably sized
    etetree = name_Internal_nodes(Tree(tree))
    
    # if internal node names would make the tree to large, dont print them 
    ascitree = etetree.get_ascii()
    ascirows = ascitree.split("\n")
    maxlen = max(map(len, ascirows))
    if maxlen > 140:
        ascitree = etetree.get_ascii(show_internal=False)

    # if the tree is still too large, limit the length of the node names
    ascirows = ascitree.split("\n")
    maxlen = max(map(len, ascirows))
    if maxlen > 140:
        for node in etetree.traverse():
            node.name = string_limit(node.name, 36)
    ascitree = etetree.get_ascii(show_internal=False)

    return ascitree

# return a list of the leaf node names from a tree
def leafname_list   (
        intree:             Tree_newick
                    ) ->    list[str]:

    tree = Tree(intree)
    leafnames = [name for name in tree.iter_leaf_names()]

    return leafnames

# return the list of the currently accepted leaf node names, without the internal populations
def accepted_leaves (
        guide_tree:         Tree_newick,
        accepted_pops:      Population_list
                    ) ->    list[str]:

    tree = Tree(guide_tree)
    tree = name_Internal_nodes(tree)
    tree.prune(accepted_pops)

    leafnames = [name for name in tree.iter_leaf_names()]

    return leafnames


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
accepted_leaf["hz_line_width"] = 6
accepted_leaf["hz_line_type"] = 1
accepted_leaf["vt_line_width"] = 2

accepted_node = NodeStyle()
accepted_node["size"] = 0
accepted_node["hz_line_width"] = 2
accepted_node["vt_line_color"] = "LimeGreen"
accepted_node["vt_line_width"] = 6
accepted_node["vt_line_type"] = 1

rejected_leaf = NodeStyle()
rejected_leaf["size"] = 0
rejected_leaf["hz_line_color"] = "Grey"
rejected_leaf["hz_line_width"] = 6
rejected_leaf["hz_line_type"] = 0
rejected_leaf["vt_line_width"] = 2

rejected_node = NodeStyle()
rejected_node["size"] = 0
rejected_node["vt_line_color"] = "Grey"
rejected_node["vt_line_width"] = 6
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

# make a tree with available distance values ultrametric, and scale the size to be easy to display
def make_ultrametric(tree):
   
    multiplier = (1/(tree.get_farthest_node()[1]))
    for node in tree.iter_descendants("postorder"):
        node.dist = node.dist*multiplier*10

    for node in tree.traverse("postorder"):
        descendants = list(node.iter_descendants("levelorder"))
        if len(descendants) > 2:
            node_1 = descendants[0]
            maxdist_1 = node_1.get_farthest_leaf()[1]
            node_2 = descendants[1]
            maxdist_2 = node_2.get_farthest_leaf()[1]
            if maxdist_1 > maxdist_2:
                node_2.dist = node_2.dist + (maxdist_1-maxdist_2)
            elif maxdist_2 > maxdist_1:
                node_1.dist = node_1.dist + (maxdist_2-maxdist_1)

    return tree

# main function implementing the drawing of the feedback tree
def visualize_decision(proposed_tree, MSC_param, BPP_outfile, proposed_changes, decision, image_name = "decision.png"):
    
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
    for node in tree.iter_descendants("preorder"):
        name = node.name
        if name in gdi_dict:
            node.add_features(gdi=f"  GDI={gdi_dict[name]}")
        if name in theta_dict:
            node.add_features(theta=f"θ={theta_dict[name]}")
        ancestor = node.up
        ancestorname = ancestor.name
        if ancestorname in tau_dict:
            node.add_features(tauval=f"tau={tau_dict[ancestorname]}")
        if ancestorname in age_dict:
            node.add_features(age=f"age={age_dict[ancestorname]}")

    for node in tree.iter_descendants("levelorder"):
        ancestor = node.up
        if ancestor.name in tau_dict:
            node.dist = tau_dict[ancestor.name]
    
    tree = make_ultrametric(tree)
    # scale theta for visualization purposes
    leaf_names = [node.name for node in tree.traverse() if node.is_leaf()]
    max_leaf_theta = 0
    for item in leaf_names:
        if item in theta_dict:
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
    age_ancestors_edited = []
    tau_ancestors_edited = []
    
    for node in tree.iter_leaves():
        name = node.name

        # write in GDI value if available
        if name in gdi_dict:
            face = AttrFace("gdi", fsize=10, fstyle="bold")
            face.hz_align = 0
            node.add_face(face, column=1, position = "aligned")

        # edit branches and circle to style and sizes requesting parameters and decision
            # edit branches and circle to accepted style
        if name in accepted_changes:
            node.set_style(accepted_leaf)
            ancestor = node.up
            ancestor.set_style(accepted_node)

            if name in theta_dict:
                nstyle = NodeStyle()
                nstyle["hz_line_width"] = 6
                nstyle["hz_line_color"] = "LimeGreen"
                nstyle["vt_line_color"] = "LimeGreen"
                nstyle["vt_line_width"] = 6
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
                nstyle["hz_line_width"] = 6
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

        # add in the tau and age values
        ancestor = node.up
        ancestorname = ancestor.name
        # if a node is the ancestor to a leaf node pair
        if len(list(ancestor.iter_descendants("levelorder"))) == 2:
            # write in tau value
            if name in left_nodes:
                if ancestorname in tau_dict and ancestorname not in tau_ancestors_edited:
                    tauface = AttrFace("tauval", fsize=10)
                    tauface.vt_align = 0
                    tauface.hz_align = 2
                    if len(age_dict) > 0:
                        tauface.margin_top = 30
                    else:
                        tauface.margin_top = 40
                    node.add_face(tauface, column=1, position = "branch-bottom")
                    tau_ancestors_edited.append(ancestorname)

            
            # write in age in generations
            if name in right_nodes:
                if ancestorname in age_dict and ancestorname not in age_ancestors_edited:
                    ageface = AttrFace("age", fsize=10, fstyle="italic")
                    ageface.vt_align = 2
                    ageface.hz_align = 0
                    ageface.margin_bottom = 20
                    node.add_face(ageface, column=1, position = "branch-top")
                    age_ancestors_edited.append(ancestorname)

        # add in actual theta value
        if name in theta_dict:
            thetaface = AttrFace("theta", fsize=10)
            thetaface.hz_align = 0
            node.add_face(thetaface, column=1, position = "branch-top")
            
    tree.render(image_name, tree_style=ts)

# graphically visualize the placement of individuals into species
def visualize_imap(current_tree, popind_dict, BPP_outfile = None, image_name = "imap.png"):
    # collect the tree
    tree = Tree(current_tree)
    if len(list(tree.iter_descendants("levelorder"))) > 1:
        tree = name_Internal_nodes(tree)
    
    # if node lengths are not available, make the tree ultrametric using the built in method
    if BPP_outfile == None:
        tree.convert_to_ultrametric()
    # if node lengths are available, make the tree ultrametric using the actual tau values
    elif BPP_outfile != None:
        # collect the split ages and add onto the tree
        tau_dict = extract_Name_TauTheta_dict(BPP_outfile)[0]
        for node in tree.iter_descendants("levelorder"):
            ancestor = node.up
            if ancestor.name in tau_dict:
                node.dist = tau_dict[ancestor.name]

        # make ultrametric
        tree = make_ultrametric(tree)

    # generate background colors
    colors = distinctipy.get_colors(len(popind_dict), pastel_factor=0.5)
    colors = [distinctipy.get_hex(color) for color in colors]
    color_dict = {pop:colors[i] for i, pop in enumerate(popind_dict)}

    # stylesheet
    nst = NodeStyle()
    nst["size"] = 0
    nst["vt_line_width"] = 2
    nst["hz_line_width"] = 2

    # add the descendand individuals, and add a colored boundary to represent the species
    for leaf in tree.traverse():
        leaf.set_style(nst)
        leafname = leaf.name
        if leafname in popind_dict:
            nst1 = NodeStyle()
            nst1["size"] = 0
            nst1["bgcolor"] = color_dict[leafname]
            nst1["vt_line_width"] = 2
            nst1["hz_line_width"] = 2
            nst1["vt_line_color"] = color_dict[leafname]
            leaf.set_style(nst1)
            individual_ids = popind_dict[leaf.name]
            for id in individual_ids:
                leaf.add_child(name = f" {id} ", dist = 0)

    # add in the name of the population to which the individuals belong to 
    edited_parents = []
    for leaf in tree.iter_leaves():
        parent = leaf.up
        if parent.name in popind_dict and parent.name not in edited_parents:
            face = TextFace(f"  {parent.name} ", fsize=14, fstyle="bold")
            face.margin_left = 10
            face.hz_align = 2
            leaf.add_face(face, column=12, position = "branch-right")
            edited_parents.append(parent.name)


    # set custom treestyle
    ts = TreeStyle()
    ts.branch_vertical_margin = 5
    ts.show_scale = False
    ts.margin_bottom = 10
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.scale = 100

    tree.render(image_name, tree_style=ts)

def visualize_progress(input_guide_tree, accepted_pops, image_name = "progress.png"):
    tree = Tree(input_guide_tree)
    tree = name_Internal_nodes(tree)
    tree.convert_to_ultrametric()
    
    # style sheet
    nstyle = NodeStyle()
    nstyle["size"] = 0
    nstyle["hz_line_color"] = "Grey"
    nstyle["vt_line_color"] = "Grey"
    nstyle["hz_line_width"] = 2
    nstyle["hz_line_type"] = 1
    nstyle["vt_line_type"] = 1
    nstyle["vt_line_width"] = 2

    nstyle2 = NodeStyle()
    nstyle2["size"] = 0
    nstyle2["hz_line_color"] = "Black"
    nstyle2["vt_line_color"] = "Black"
    nstyle2["hz_line_width"] = 4
    nstyle2["hz_line_type"] = 0
    nstyle2["vt_line_type"] = 0
    nstyle2["vt_line_width"] = 4

    nstyle3 = NodeStyle()
    nstyle3["size"] = 0
    nstyle3["hz_line_color"] = "Black"
    nstyle3["vt_line_color"] = "Grey"
    nstyle3["hz_line_width"] = 4
    nstyle3["hz_line_type"] = 0
    nstyle3["vt_line_type"] = 1
    nstyle3["vt_line_width"] = 2

    # set basic style
    for node in tree.traverse():
        node.set_style(nstyle)

    # set styles according to acceptance status, and add extra annotation for node names
    for node in tree.traverse():
        if node.name in accepted_pops:
            node.set_style(nstyle2)
            face = TextFace(f"{string_limit(node.name, 8)}", fsize=14, fstyle="bold")
            face.hz_align = 2
            face.margin_right = 20
            node.add_face(face, column=1, position = "branch-top")
        elif node.is_leaf() == False:
            face = TextFace(f"{string_limit(node.name, 8)}", fsize=8)
            face.hz_align = 2
            face.margin_right = 20
            node.add_face(face, column=1, position = "branch-bottom")

    for node in tree.traverse():
        if node.name not in accepted_pops:
            parent = node.up
            if parent.name in accepted_pops:
                parent.set_style(nstyle3)

    tree.render(image_name, tree_style=ts)

# visualize the output produced in the starting topology section
def visualize_tree(input_tree, image_name = "tree.png"):
    tree = Tree(input_tree)
    tree.convert_to_ultrametric()
    
    # style sheet
    nstyle = NodeStyle()
    nstyle["size"] = 0
    nstyle["hz_line_color"] = "Black"
    nstyle["vt_line_color"] = "Black"
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2
    nstyle["hz_line_type"] = 0
    nstyle["vt_line_type"] = 0
    
    # set basic style
    for node in tree.traverse():
        node.set_style(nstyle)

    tree.render(image_name, tree_style=ts)