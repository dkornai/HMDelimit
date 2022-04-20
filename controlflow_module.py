from helper_functions import read_MasterControl
from bpp_cfile_module import get_known_BPP_param

# find trees in all possible locations
def get_Tree(input_mc_file):
    mc_dict = read_MasterControl(input_mc_file)

    trees = {"tree_HM"        :mc_dict["tree_HM"], 
             "tree_start"     :mc_dict["tree_start"], 
             "tree_BPPA00"    :get_known_BPP_param(mc_dict, "A00" )["newick"],
             "tree_BPPA11"    :get_known_BPP_param(mc_dict, "A11" )["newick"],
             "tree_BPPA01"    :get_known_BPP_param(mc_dict, "A01" )["newick"],}

    return trees

# find imap files in all possible locations
def get_Imap(input_mc_file):
    mc_dict = read_MasterControl(input_mc_file)

    imapfiles = {"imap_MCF"   :mc_dict["file_imap"],
                 "imap_hm"    :get_known_BPP_param(mc_dict, "A00" )["Imapfile"],
                 "imap_delim" :get_known_BPP_param(mc_dict, "A11" )["Imapfile"],
                 "imap_start" :get_known_BPP_param(mc_dict, "A01" )["Imapfile"],}

    return imapfiles

# find alignment files in all possible locations
def get_Imap(input_mc_file):
    mc_dict = read_MasterControl(input_mc_file)

    seqfiles = {"seq_MCF"     :mc_dict["file_seq"],
                "seq_hm"      :get_known_BPP_param(mc_dict, "A00" )["seqfile"],
                "seq_delim"   :get_known_BPP_param(mc_dict, "A11" )["seqfile"],
                "seq_start"   :get_known_BPP_param(mc_dict, "A01" )["seqfile"],}

    return seqfiles


# scan the tree presence dict to decide which state the program will start from
def find_initial_State(input_mc_file):
    print("\nCONTROL FLOW:\n")

    mc_dict = read_MasterControl(input_mc_file)
    tree_state = get_Tree(input_mc_file)
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

    print("\t\t-- PROCEEDING TO NEXT STEP --")

    return p_state


