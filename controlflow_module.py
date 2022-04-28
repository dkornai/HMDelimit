'''
THIS MODULE DECIDES WHICH STEPS THE PIPELINE SHOULD TAKE, 
DEPENDING ON THE DATA SUPPLIED BY THE USER. THE MODULE
ALSO PERFORMS THE APPROPRIATE TYPES OF CHECKS ON DATA AND
PARAMETER VALIDITY
'''
## DEPENDENCIES
# CHECK HELPER FUNCTIONS
from check_helper_functions import check_BPP_ctl_validity

# CHECK PARAM MODULE
from check_param_module import check_BPP_param
from check_param_module import check_A11_input_compat
from check_param_module import check_A00_input_compat

# HELPER FUNCTIONS
from helper_functions import read_MasterControl

# BPP CFILE MODULE
from bpp_cfile_module import get_known_BPP_param
from bpp_cfile_module import get_user_BPP_param

## DATA DEPENDENCIES
from data_dicts import clprnt

## TYPE HINTS
from custom_types import Master_control_dict
from custom_types import Master_control_file
from custom_types import Tree_newick

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

# scan the master control file and BPP control file to determine which stages to execute
def find_initial_State  (
        input_mc_file:          Master_control_file
                        ):

    print(f"\n{clprnt.BLUE}CHECKING OF CONTROL FLOW{clprnt.end}\n")

    mc_dict = read_MasterControl(input_mc_file)
    tree_state = get_Tree(mc_dict)
    
    ## DECIDE WHICH STATE TO ENTER INTO, DEPENDING ON THE PARAMETERS THAT WERE PROVIDED

    # if a guide tree for the HM stage is specified in it's control file, skip all stages before HM
    if tree_state["tree_BPPA00"] != "?":
        p_state = "A00"
        print(f"HM guide tree identified in the BPP control file: {mc_dict['ctl_file_HM']}\n\n\tTree: {tree_state['tree_BPPA00']}")
    
    # if a guide tree for the HM stage is specified in the MCF, skip all stages before HM
    elif tree_state["tree_HM"] != "?":
        p_state = "A00"
        print(f"HM guide tree identified in Master Control File!\n\n\tTree: {tree_state['tree_HM']}")

    # if a tree for the SD stage is specified in a seperate control file
    elif tree_state["tree_BPPA11"] != "?":
        p_state = "A11+A00"
        print(f"BPP A11 guide tree identified in the BPP control file: {mc_dict['ctl_file_delim']}\n\n\tTree: {tree_state['tree_BPPA11']}")
    # if a starting tree for A11 is specified in the MCF
    elif tree_state["tree_start"] != "?":
        p_state = "A11+A00"
        print(f"BPP A11 guide tree identified in Master Control File!\n\n\tTree: {tree_state['tree_start']}")

    # if a guide tree for the starting topolgy inference stage is specified the control file, and the A11 option is True
    elif tree_state["tree_BPPA01"] != "?" and mc_dict["execute_A11"] == "True":
        p_state = "A01+A11+A00"
        print(f"Starting Topology inference helper tree identified in the BPP control file: {mc_dict['ctl_file_delim']}\n\n\tTree: {tree_state['tree_BPPA01']}")
        print("Execute A11 option set to 'True' in Master Control File")

    # if a guide tree for the starting topolgy inference stage is specified the control file, and an A11 control file is specified
    elif tree_state["tree_BPPA01"] != "?" and mc_dict["ctl_file_delim"] != "?":
        p_state = "A01+A11+A00"
        print(f"Starting Topology inference helper tree identified in the BPP control file: {mc_dict['ctl_file_delim']}\n\n\tTree: {tree_state['tree_BPPA01']}")
        print(f"BPP A11 control file specified in Master Control File: {mc_dict['ctl_file_delim']}")
    
    # if a guide tree for the starting topolgy inference stage is specified the control file,
    elif tree_state["tree_BPPA01"] != "?":
        p_state = "A01+A00"
        print(f"Starting Topology inference helper tree identified in the BPP control file: {mc_dict['ctl_file_delim']}\n\n\tTree: {tree_state['tree_BPPA01']}")

    # if no tree is found, but the user wants the A11 step as indicated by the execute_A11 option
    elif mc_dict["execute_A11"] == "True":
        p_state = "A01+A11+A00"
        print("No trees found in the Master Control File, or associated BPP control file")
        print("The execute A11 option is set to 'True'.")

    # if no tree is found, but the user wants the A11 step as indicated by a specified control file
    elif mc_dict["ctl_file_delim"] != "?":
        p_state = "A01+A11+A00"
        print("No trees found in the Master Control File, or associated BPP control file")
        print(f"BPP A11 control file specified in Master Control File: {mc_dict['ctl_file_delim']}")
    
    # if no tree is found
    else:
        p_state = "A01+A00"
        print("No trees found in the Master Control File, or associated BPP control files.")

    ## PROVIDE USER FEEDBACK
    print("\nAccordingly, the program will:\n")
    if p_state == "A00":
        print("1) Use BPP A00 and the the Hierarchical Method to find the optimal species delimitation\ngiven the guide tree and imap provided by the user")
    
    if p_state == "A01+A00":
        print("1) Use BPP A01 to infer a species phylogeny from imap and seqfile")
        print("2) Use BPP A00 and the the Hierarchical Method to find the optimal species delimitation\ngiven the guide tree produced in the previous A01 stage")

    if p_state == "A11+A00":
        print("1) Use BPP A11 to infer a species delimitation from the starting tree and imap provided by the user")
        print("2) Use BPP A00 and the the Hierarchical Method to further refine the delimitation produced by the previous A11 stage")
    
    if p_state == "A01+A11+A00":
        print("1) Use BPP A01 to infer a species phylogeny from imap and seqfile")
        print("1) Use BPP A11 to infer a species delimitation given the starting tree produced in the previous A01 stage")
        print("3) Use BPP A00 and the the Hierarchical Method to further refine the delimitation produced by the previous A11 stage")

    print()

    return p_state

# check if the BPP control files that are supplied by the user, and relevant to the execution only contain valid paramters
def controlled_BPP_cfile_check  (
        input_mc_file:                  Master_control_file, 
        p_state:                        str
                                ):

    mc_dict = read_MasterControl(input_mc_file)
    
    A01_ctl = mc_dict['ctl_file_phylo']
    A11_ctl = mc_dict['ctl_file_delim']
    A00_ctl = mc_dict['ctl_file_HM']
    
    ## P STATE DEPENDENT CHECKING
    all_compatible = []
        
    if   p_state == "A00":
        if A00_ctl != "?":
            print(f"{clprnt.BLUE}INITAL CHECK OF USER SUPPLIED BPP CONTROL FILES{clprnt.end}")
        
        if A00_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A00 CONTROL FILE: {A00_ctl}")
            if check_BPP_ctl_validity(A00_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)

    elif p_state == "A01+A00":
        if A11_ctl != "?" or A00_ctl != "?":
            print(f"{clprnt.BLUE}INITAL CHECK OF USER SUPPLIED BPP CONTROL FILES{clprnt.end}")
        
        if A01_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A11 CONTROL FILE: {A01_ctl}")
            if check_BPP_ctl_validity(A01_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)
        if A00_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A00 CONTROL FILE: {A00_ctl}")
            if check_BPP_ctl_validity(A00_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)

    elif p_state == "A11+A00":
        if A11_ctl != "?" or A00_ctl != "?":
            print(f"{clprnt.BLUE}INITAL CHECK OF USER SUPPLIED BPP CONTROL FILES{clprnt.end}")
        
        if A11_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A11 CONTROL FILE: {A11_ctl}")
            if check_BPP_ctl_validity(A11_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)
        if A00_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A00 CONTROL FILE: {A00_ctl}")
            if check_BPP_ctl_validity(A00_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)
        
    elif p_state == "A01+A11+A00":
        if A11_ctl != "?" or A00_ctl != "?" or  A01_ctl != "?":
            print(f"{clprnt.BLUE}INITAL CHECK OF USER SUPPLIED BPP CONTROL FILES{clprnt.end}")
        
        if A01_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A01 CONTROL FILE: {A01_ctl}")
            if check_BPP_ctl_validity(A01_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)
        if A11_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A11 CONTROL FILE: {A11_ctl}")
            if check_BPP_ctl_validity(A11_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)
        if A00_ctl != "?":
            print(f"\n\tCHECKING SUPPLIED A00 CONTROL FILE: {A00_ctl}")
            if check_BPP_ctl_validity(A00_ctl):
                all_compatible.append(True)
            else:
                all_compatible.append(False)

    ## FINAL DECISION
    if   all(comp == True for comp in all_compatible) and len(all_compatible) > 0:
        print("\n\t[*] All supplied BPP control files are valid")
    
    elif any(comp == False for comp in all_compatible):
        print("[X] ERROR: BPP CONTROL FILE(S) CONTAIN UNKNOWN PARAMETERS!")
        exit()


# perform the appropriate checks for BPP parameters in the stages that will be executed
def controlled_BPP_parameter_check  (
        input_mc_file:                      Master_control_file, 
        p_state:                            int
                                    ):

    mc_dict = read_MasterControl(input_mc_file)
    compatible = False

    ## PERFORM MODE DEPENDENT CHECKS
    if   p_state == "A00":
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00")
        
        if check_BPP_param(A00_param, A00_source, "A00"):
            compatible = True

    elif p_state == "A01+A00":
        A01_param, A01_source = get_user_BPP_param(mc_dict, "A01")
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00")
        
        if check_BPP_param(A01_param, A01_source, "A01"):
            if check_A00_input_compat(A01_param, A00_param):
                if check_BPP_param(A00_param, A00_source, "A00"):
                    compatible = True

    elif p_state == "A11+A00":
        A11_param, A11_source = get_user_BPP_param(mc_dict, "A11")
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00", after_A11 = True)
        
        if check_BPP_param(A11_param, A11_source, "A11"):
            if check_A00_input_compat(A11_param, A00_param):
                if check_BPP_param(A00_param, A00_source, "A00", after_A11 = True):
                    compatible = True

    elif p_state == "A01+A11+A00":
        A01_param, A01_source = get_user_BPP_param(mc_dict, "A01")
        A11_param, A11_source = get_user_BPP_param(mc_dict, "A11")
        A00_param, A00_source = get_user_BPP_param(mc_dict, "A00", after_A11 = True)
        
        if check_BPP_param(A01_param, A01_source, "A01"):
            if check_A11_input_compat(A01_param, A11_param):
                if check_BPP_param(A11_param, A11_source, "A11"):
                    if check_A00_input_compat(A11_param, A00_param):
                        if check_BPP_param(A00_param, A00_source, "A00", after_A11 = True):
                            compatible = True

    ## EXIT PROGRAM IF INCOMPATIBILITIES ARE FOUND
    if compatible == False:
        exit()