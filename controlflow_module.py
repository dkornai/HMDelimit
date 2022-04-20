from check_param_module import check_A01_to_A11_compatibility, check_A11_to_A00_compatibility, check_BPP_cfile
from custom_types import Master_control_dict, Tree_newick
from helper_functions import read_MasterControl
from bpp_cfile_module import get_known_BPP_param
from bpp_cfile_module import get_user_BPP_param

from data_dicts import col_print

# find trees in all possible locations
def get_Tree(
        mc_dict:    Master_control_dict
            ) ->    dict[str, Tree_newick]:

    trees = {"tree_HM"        :mc_dict["tree_HM"], 
             "tree_start"     :mc_dict["tree_start"], 
             "tree_BPPA00"    :get_known_BPP_param(mc_dict, "A00" )["newick"],
             "tree_BPPA11"    :get_known_BPP_param(mc_dict, "A11" )["newick"],
             "tree_BPPA01"    :get_known_BPP_param(mc_dict, "A01" )["newick"],}

    return trees

# scan the tree presence dict to decide which state the program will start from
def find_initial_State(input_mc_file):
    print(f"\n{col_print.BLUE}CHECKING OF CONTROL FLOW{col_print.RESETC}\n")

    mc_dict = read_MasterControl(input_mc_file)
    tree_state = get_Tree(mc_dict)
    
    # program states are entered in this order, depending on whether a tree for a given stage was specified
        # if a guide tree for the HM stage is specified, skip all stages before HM
    if tree_state["tree_HM"] != "?":
        p_state = 3
        print(f"HM guide tree identified in Master Control File!\n\n\tTree: {tree_state['tree_HM']}")
        # if a starting tree for A11 is specified, skip all other stages before
    elif tree_state["tree_start"] != "?":
        p_state = 2
        print(f"Starting Delimitation guide tree identified in Master Control File!\n\n\tTree: {tree_state['tree_start']}")
        # if a guide tree for the HM stage is specified in it's control file, and no trees are in the MCF
    elif tree_state["tree_BPPA00"] != "?":
        p_state = 3
        print(f"HM guide tree identified in the BPP control file: {mc_dict['ctl_file_HM']}\n\n\tTree: {tree_state['tree_BPPA00']}")
        # if a tree for the SD stage is specified in it's control file, and no trees are in the MCF
    elif tree_state["tree_BPPA11"] != "?":
        p_state = 2
        print(f"Starting delimitation guide tree identified in the BPP control file: {mc_dict['ctl_file_delim']}\n\n\tTree: {tree_state['tree_BPPA11']}")
        # if a guide tree for the starting topolgy inference stage is specified in it's control file, and no trees are in the MCF
    elif tree_state["tree_BPPA01"] != "?":
        p_state = 1
        print(f"Starting Topology inference helper tree identified in the BPP control file: {mc_dict['ctl_file_delim']}\n\n\tTree: {tree_state['tree_BPPA01']}")
        # if no tree is found
    else:
        p_state = 1
        print("No trees found in the Master Control File, or associated BPP control files.")

    print("\nAccordingly, the program will:\n")

    # print user feedback about what the program is going to do
    stage_1_desc = "Use BPP A01 to infer a starting phylogenetic tree\n"
    stage_2_desc = "Use BPP A11 to infer a starting delimitation and an associated guide tree\n"
    stage_3_desc = "Use BPP A00 and the the Hierarchical Method to find the optimal species delimitation\n   given the starting delimitation and the guide tree\n"
    
    if p_state == 1:
        print("1)", stage_1_desc)
        print("2)", stage_2_desc)
        print("3)", stage_3_desc)

    if p_state == 2:
        print("1)", stage_2_desc)
        print("2)", stage_3_desc)
    
    if p_state == 3:
        print("1)", stage_3_desc)

    return p_state

# perform the appropriate checks for BPP parameters in the stages that will be executed
def controlled_check(input_mc_file, p_state):

    mc_dict = read_MasterControl(input_mc_file)
    compatible = False

    if   p_state == 3:
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00")
        
        if check_BPP_cfile(A00_param, A00_source, "A00"):
            compatible = True

    elif p_state == 2:
        A11_param, A11_source = get_user_BPP_param(mc_dict, "A11")
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00", mask_A11 = True)
        
        if check_BPP_cfile(A11_param, A11_source, "A11"):
            if check_A11_to_A00_compatibility(A11_param, A00_param):
                if check_BPP_cfile(A00_param, A00_source, "A00"):
                    compatible = True

    elif p_state == 1:
        A01_param, A01_source = get_user_BPP_param(mc_dict, "A01")
        A11_param, A11_source = get_user_BPP_param(mc_dict, "A11")
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00", mask_A11 = True)
        
        if check_BPP_cfile(A01_param, A01_source, "A01"):
            if check_A01_to_A11_compatibility(A01_param, A11_param):
                if check_BPP_cfile(A11_param, A11_source, "A11"):
                    if check_A11_to_A00_compatibility(A11_param, A00_param):
                        if check_BPP_cfile(A00_param, A00_source, "A00"):
                            compatible = True

    return compatible