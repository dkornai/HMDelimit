## DEPENDENCDIES

# PYTHON STANDARD DEPENDENCIES
import re
import copy
import warnings

# PYTHON LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

# HELPER FUNCTION DEPENDENCIES
from helper_functions import Imap_to_PopInd_Dict
from helper_functions import Imap_to_List

# TYPE HINTING DEPENDENCIES
from custom_types import BPPControl_dict
from custom_types import Imap_list
from custom_types import Population_list
from custom_types import Newick_tree
from custom_types import Imap_file

## ENCODING FUNCTIONS

# generate unique IDs that are used to avoid possible overlaps between population names in the user supplied Imap
'''
This is necessary as BPP A11 does not return the imap corresponding to the final delimitation, 
only the newick tree, and the list of accepted populations. This is a problem, because A11 
might also change the tree topology, and not just merge populations according to the original 
tree. This becomes a problem when we want to generate the new imap, which includes the population
merges implemented in A11.
'''
def generate_nonoverlap_IDs(imapfile:Imap_file):
    pops = list(Imap_to_PopInd_Dict(imapfile).keys())
    real_tempname_dict = {}
    for i, pops in enumerate(pops):
        # replace original name with a nonoverlap numeric ID that is always formatted identically
        real_tempname_dict[pops] = f"PN_{str(i).zfill(3)}"

    return real_tempname_dict

# remap the population names in the BPP control file to their unqie IDs
def enforce_no_overlaps_BPP_cfile(input_BPP_cdict:BPPControl_dict, real_tempname_dict:dict):
    cdict = copy.deepcopy(input_BPP_cdict)

    # The names on the s&t line can be subbed using regex, as there is always a whitespace between the names
    s_and_t = cdict["species&tree"]
        # split off the number of species, as this may confuse the regex if 
        # any of the species names has an identical number
    s_and_t = f"{re.split('[0-9]+', s_and_t, maxsplit = 1)[1]}"
    n = re.sub(f"{s_and_t}","", cdict["species&tree"])
    s_and_t += "  " # pad with extra spaces to avoid bugs where the last species is not subbed
    
    for realname in real_tempname_dict:
        s_and_t = re.sub(f"\s{realname}\s", f" {real_tempname_dict[realname]} ", s_and_t)
    
    cdict["species&tree"] = f"{n}{s_and_t}"
    
    # The names in the tree need to be subbed in a tree data structure, as the lack of whitespaces
    # in the classic newick format means that names can be hard to separate 
    newick = cdict["newick"]
    tree = Tree(newick)
    for node in tree.traverse("postorder"):
        for realname in real_tempname_dict:
            if node.name == realname:
                node.name = real_tempname_dict[realname]
                break
    
    cdict["newick"] = tree.write(format=9)

    return cdict


# remap the population names in the imap file to their unique IDs
def enforce_no_overlaps_Imap(imapfile:Imap_file, real_tempname_dict:dict):
    imaplist = Imap_to_List(imapfile)
    
    # scan through the populations column and replace population names with their non-overlapping IDs
    output = copy.deepcopy(imaplist)
    for i, popname in enumerate(imaplist[1]):
        for realname in real_tempname_dict:
            if popname == realname:
                output[1][i] = real_tempname_dict[popname]
                break
    
    return output


# final wrapper function intended to implement unique ID encoding for the BPP control file and Imap simultaneously
def uniqueID_encoding(input_BPP_cdict:BPPControl_dict, imap_name:str):
    remap_dict = generate_nonoverlap_IDs(input_BPP_cdict["Imapfile"])
    BPP_cdict  = enforce_no_overlaps_BPP_cfile(input_BPP_cdict, remap_dict)
    imaplist   = enforce_no_overlaps_Imap(input_BPP_cdict["Imapfile"], remap_dict)
    BPP_cdict["Imapfile"] = imap_name

    return BPP_cdict, imaplist, remap_dict


## DECODING FUNCTIONS

# do a remap on merged population IDs using a dict that only includes population ID components
def componentwise_remap(label, remap_dict:dict):
    # split into 6 character chunks, to avoid possible contamination by overlaps
    label_components = [label[i:i+6] for i in range(0, len(label), 6)]
    # replace using the dict, component by component
    label_components = [component.replace(component, remap_dict[component]) for component in label_components]
    
    return "".join(label_components)


# get the output of BPP A11 (which had the unique IDs injected), and infer the Imap that would
# be produced using the original naming scheme of the user.
def infer_Imap_from_A11_output(imap_list:Imap_list, new_pops:Population_list, real_tempname_dict:dict):
    imaplist = imap_list
    unique_pops = imaplist[1]

    # generate mapping of accepted (possibly) merged populations to old single populations
    merged_to_single_dict = {}
    for old in unique_pops:
        for new in new_pops:
            if old in new:
                merged_to_single_dict[old] = new

    # remap the single populations to their new (possibly) merged versions, which are listed in "new_pops"
    new_pops = [merged_to_single_dict[popid] for popid in unique_pops]
    imaplist[1] = new_pops

    # remap the unique IDs to their original user supplied names
        # invert the remapping dict
    tempname_real_dict = dict((v, k) for k, v in real_tempname_dict.items())    
        # remap the IDs using the dict
    for x, label in enumerate(imaplist[1]):
        imaplist[1][x] = componentwise_remap(label, tempname_real_dict)
               
    return imaplist


# remap the unique IDs injected into A11 to their original names, this is challenging due to merged names
def infer_Tree_from_A11_output(newick:Newick_tree, real_tempname_dict:dict) -> Newick_tree:
    tempname_real_dict = dict((v, k) for k, v in real_tempname_dict.items())
    
    tree = Tree(newick)
    for node in tree.traverse("postorder"):
        if node.name != "?":
            node.name = componentwise_remap(node.name, tempname_real_dict)
    
    return tree.write(format=9)


# final wrapper function to decode the machine generated unique IDs into the original user supplied ones
def uniqueID_decoding(imap_list:Imap_list, new_pops:Population_list, newick_tree:Newick_tree, remap_dict:dict):
    imap = infer_Imap_from_A11_output(imap_list, new_pops, remap_dict)
    tree = infer_Tree_from_A11_output(newick_tree, remap_dict)

    return tree, imap


### TESTS 
# generate_nonoverlap_IDs("D_L10_imap.txt")
# generate_nonoverlap_IDs("D_TMS_imap.txt")
# generate_nonoverlap_IDs("D_TMS_imap_mod.txt")
# generate_nonoverlap_IDs("D_ROT_imap.txt")

# test_cfile_dict = {"species&tree": "66  D  C6 B BedBug Bug EEEBug F",
#                   "newick"      :"((((D,C6),(B,BedBug)),Bug), (EEEBug, F));"}     
 
# print(enforce_no_overlaps_BPP_cfile(test_cfile_dict, generate_nonoverlap_IDs("D_TMS_imap_mod.txt")))

# print(enforce_no_overlaps_Imap("D_TMS_imap_mod.txt", generate_nonoverlap_IDs("D_TMS_imap_mod.txt")))

# print(infer_Imap_from_A11_output(Imap_to_List("D_TMS_imap_mod_un.txt"), 
#                                  ['PN_001', 'PN_002', 'PN_003','PN_004', 'PN_000', 'PN_005', 'PN_006'], 
#                                  generate_nonoverlap_IDs("D_TMS_imap_mod.txt")))
# print(infer_Imap_from_A11_output(Imap_to_List("D_TMS_imap_mod_un.txt"), 
#                                  ['PN_001', 'PN_002', 'PN_003PN_004', 'PN_000', 'PN_005PN_006'], 
#                                  generate_nonoverlap_IDs("D_TMS_imap_mod.txt")))
