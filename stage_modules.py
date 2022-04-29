## DEPENDENCIES
# STANDARD LIBRARY
import os
import shutil

# HELPER FUNCTIONS
from helper_functions import Imap_to_IndPop_Dict 
from helper_functions import Imap_to_PopInd_Dict
from helper_functions import create_TargetDir
from helper_functions import pretty
from helper_functions import Imap_to_List
from helper_functions import get_HM_parameters
from helper_functions import dict_to_bppcfile
from helper_functions import list_To_Imap
from helper_functions import read_MasterControl
from helper_functions import BPP_run
from helper_functions import extract_Speciestree
from helper_functions import extract_Pops
from helper_functions import string_limit

# TREE HELPER FUNCTIONS
from tree_helper_functions import visualize_progress, write_Tree

# BPP CONTROL FILE RELATED FUNCTIONS
from bpp_cfile_module import get_known_BPP_param 
from bpp_cfile_module import generate_unkown_BPP_param
from bpp_cfile_module import generate_unknown_BPP_tree
from bpp_cfile_module import proposal_compliant_BPP_param

# UNIQUE ID ENCODING AND DECONDING FUNCTIONS
from uniqueID_module import uniqueID_encoding
from uniqueID_module import uniqueID_decoding 

# PROPOSAL FUNCTIONS
from proposal_module import HMproposal
from proposal_module import get_HM_StartingState
from proposal_module import get_HM_results

# DECISION FUNCTIONS
from decision_module import decisionModule
from decision_module import get_MSC_param

# DEEP CONFLICT CHECKING FUNCTIONS
from check_conflict_functions import check_GuideTree_Imap_MSA_compat

# TREE HELPER FUNCTIONS
from tree_helper_functions import tree_ASCII
from tree_helper_functions import visualize_decision
from tree_helper_functions import visualize_imap
from tree_helper_functions import leafname_list
from tree_helper_functions import accepted_leaves

## DATA DEPENDENCIES
from data_dicts import clprnt

## TYPE HINTS
from custom_types import Imap_list
from custom_types import Tree_newick
from custom_types import Population_list
from custom_types import Master_control_file


# infer the starting topology before any delimitation steps
def StartingTopolgy (
        input_mcfile:       Master_control_file
                    ) ->    Tree_newick:

    parent_dir = os.getcwd()
    print(f"{clprnt.BLUE}\nBEGINNING STARTING PHYLOGENY INFERENCE\n{clprnt.end}")
    
    # get master control file parameters
    mc_dict = read_MasterControl(input_mcfile)
    
    # create the target directory specific to the step, and the name of the MCF
    target_dir = f'{input_mcfile[0:-4]}_0_StartPhylo'
    create_TargetDir(target_dir)

    # set up the BPP control file specific to the A01 stage
    BPP_A01_cfile_name = "BPP_A01_StartPhylo.ctl"
    BPP_cdict = get_known_BPP_param(input_mc_dict = mc_dict, BPP_mode = 'A01')
    BPP_cdict = generate_unkown_BPP_param(BPP_cdict)
    BPP_cdict = generate_unknown_BPP_tree(BPP_cdict)
    # print feedback to the user
    print("\nBPP CONTROL FILE:")
    pretty(BPP_cdict)

    # write the relevant files to the target directory
    dict_to_bppcfile    (BPP_cdict, os.path.join(target_dir, BPP_A01_cfile_name))
    shutil.copy         (src = BPP_cdict['seqfile'],  dst = target_dir)
    shutil.copy         (src = BPP_cdict['Imapfile'], dst = target_dir)
    
        # STARTING PHYLOGENY #
    ###############################
    #-----------------------------#
    os.chdir(target_dir)
    
    # run bpp
    BPP_run(BPP_A01_cfile_name)
        
    # extract the species tree
    tree = extract_Speciestree(BPP_A01_cfile_name)

    # write resulting tree for manual inspection
    write_Tree(tree, "OUTPUT_TREE.txt")

    os.chdir(parent_dir)
    #-----------------------------#
    print("\n>> RESULTS OF STARTING PHYLOGENY INFERENCE:")
    print(f"\t\t\nSTARTING NEWICK TREE:\n\n\t{string_limit(tree, 96)}")
    print(f"\nSTARTING ASCII TREE:")
    print(tree_ASCII(tree))

    return tree


# infer the starting delimitation. This consists of a guide tree and an associated Imap
def StartingDelimitation(
        input_mcfile:           Master_control_file, 
        starting_tree:          Tree_newick = None
                        ) ->    tuple[Tree_newick, Imap_list]:

    parent_dir = os.getcwd()
    print(f"{clprnt.BLUE}\nBEGINNING STARTING DELIMITATION{clprnt.end}\n")

    # get master control file parameters
    mc_dict = read_MasterControl(input_mcfile)
    
    # create the target directory specific to the step
    target_dir = f'{input_mcfile[0:-4]}_1_StartDelim'
    create_TargetDir(target_dir)

    # set up the BPP control file specific to the A11 stage
    BPP_A11_cfile_name = "BPP_A11_StartDelim.ctl"
    BPP_cdict = get_known_BPP_param(input_mc_dict = mc_dict, BPP_mode = 'A11')
    BPP_cdict = generate_unkown_BPP_param (BPP_cdict)
    # overwrite any existing starting tree if one was generated in the A01 step or supplied in the MCF
    if starting_tree != None:
        BPP_cdict["newick"] = starting_tree
    # print feedback to the user
    print("\nBPP CONTROL FILE:")
    pretty(BPP_cdict)

    # unique ID encoding
    imap_unique_ids_name = "Imap_UniqueID.txt"
    BPP_cdict, imap_unique_ids, remap_dict = uniqueID_encoding(BPP_cdict, imap_unique_ids_name)

    # write the relevant files to the target directory
    dict_to_bppcfile    (BPP_cdict, os.path.join(target_dir, BPP_A11_cfile_name))
    list_To_Imap        (imap_unique_ids, os.path.join(target_dir, imap_unique_ids_name))
    shutil.copy         (src = BPP_cdict['seqfile'],  dst = target_dir)

       # STARTING DELIMITATION #
    ###############################
    #-----------------------------#
    os.chdir(target_dir)
    
    # run BPP
    BPP_run(BPP_A11_cfile_name)
        
    # capture output (encoded with unique IDs)
    tree_encoded = extract_Speciestree(BPP_A11_cfile_name)
    pops_encoded = extract_Pops(BPP_A11_cfile_name)
        
    # decode unique IDs back to user supplied names
    guide_tree, imap = uniqueID_decoding(imap_unique_ids, pops_encoded, tree_encoded, remap_dict)
        
    # write resulting Imap and tree for manual inspection if needed
    list_To_Imap(imap, "OUTPUT_IMAP.txt")
    write_Tree(guide_tree, "OUTPUT_TREE.txt")

    os.chdir(parent_dir)
    #-----------------------------#
    ###############################

    print("\n>> RESULTS OF THE STARTING DELIMITATION STAGE:")
    print(f"\t\t\nGUIDE TREE:\n\n\t{string_limit(guide_tree, 96)}")
    print(f"\nGUIDE ASCII TREE:")
    print(tree_ASCII(guide_tree))
    print(f"\t\t\nIMAP:\n")
    pretty(Imap_to_PopInd_Dict(imap))

    return guide_tree, imap, 


# perform one iteration of the hierarchical method.
def HMIteration (
        input_mcfile:           Master_control_file, 
        input_guide_tree:       Tree_newick, 
        input_indpop_dict, 
        input_accepted_pops:    Population_list, 
        halt_pop_number:        int, 
        step:                   int
                ) ->            tuple[Population_list, bool]:

    parent_dir = os.getcwd()
    print(f"{clprnt.BLUE}\nBEGINNING ITERATION {step} OF THE HIERARCHICAL METHOD{clprnt.end}\n")
    
    # get master control file parameters
    mc_dict = read_MasterControl(input_mcfile)

    # get HM decision specific parameters (eg GDI threshreal, mutation rate...)
    hm_param = get_HM_parameters(input_mc_dict = mc_dict)

    # create the target directory specific to the step
    target_dir = f'{input_mcfile[0:-4]}_2_HM_{step}'
    create_TargetDir(target_dir_name=target_dir)
    
    # generate a proposal based on the previously accepted results
    prop_change, prop_tree, prop_imap = HMproposal(guide_tree_newick = input_guide_tree,
                                                   base_indpop_dict  = input_indpop_dict,
                                                   current_pops_list = input_accepted_pops,
                                                   mode              = hm_param["mode"])
    prop_imap_name = "proposed_imap.txt"
    
    # set up the control file specific to the A00 stage
    proposed_cfile_name = f"BPP_A00_HM_{step}.ctl"
    BPP_cdict = get_known_BPP_param(input_mc_dict = mc_dict, BPP_mode = 'A00')
    BPP_cdict = generate_unkown_BPP_param(BPP_cdict) 
    BPP_cdict = proposal_compliant_BPP_param(BPP_cdict, prop_imap, prop_imap_name, prop_tree)
    # print feedback to the user
    print("\nBPP CONTROL FILE:\n")
    pretty(BPP_cdict)
    
    # write the relevant files
    list_To_Imap        (prop_imap, os.path.join(target_dir, prop_imap_name))
    shutil.copy         (src = BPP_cdict['seqfile'], dst = target_dir)
    dict_to_bppcfile    (BPP_cdict, os.path.join(target_dir, proposed_cfile_name))

           # HM ITERATION #
    ###############################
    #-----------------------------#
    os.chdir(target_dir)
        
    # run BPP 
    BPP_run(proposed_cfile_name)
    
    # make decision about which proposals to accept based on BPP results and HM decision criteria
    accepted, to_iterate, decision = decisionModule(hm_param         = hm_param,
                                                    BPP_outfile      = os.path.join(BPP_cdict["outfile"]),
                                                    proposed_changes = prop_change,
                                                    accepted_pops    = input_accepted_pops,
                                                    halt_pop_number  = halt_pop_number)

    # write tree and imap corresponding to results
    imap, tree = get_HM_results(input_guide_tree, input_indpop_dict, accepted)
    list_To_Imap    (imap, "OUTPUT_IMAP.txt")
    write_Tree      (tree, "OUTPUT_TREE.txt")
    
    # write images showing how the delimitation has been impacted by the current step
    visualize_decision  (prop_tree, get_MSC_param(BPP_cdict["outfile"], prop_change, hm_param), BPP_cdict["outfile"], prop_change, decision)
    visualize_imap      (tree, Imap_to_PopInd_Dict(imap), BPP_cdict["outfile"])
    visualize_progress  (input_guide_tree, accepted)

    os.chdir(parent_dir)
    #-----------------------------#
    ###############################

    print(f"\n>> RESULTS AFTER ITERATION {step}:\n")
    print("THE CURRENTLY ACCEPTED SPECIES ARE:\n")
    for population in leafname_list(tree):
        print(f"\t{population}")
    print()
    print(f"CURRENT TREE:\n\n\t{string_limit(tree, 96)}")
    print(f"\nCURRENT ASCII TREE:")
    print(tree_ASCII(tree))
    print(f"\nCURRENT IMAP:\n")
    pretty(Imap_to_PopInd_Dict(imap))

    return accepted, to_iterate


# final wrapper function for starting and iterating through the Hierarchical Method
def HierarchicalMethod  (
        input_mcfile:       Master_control_file, 
        input_guide_tree:   Tree_newick = None, 
        input_imap:         Imap_list = None,
                        ):
    
    print(f"\n{clprnt.BLUE}BEGINNING THE HIERARCHICAL METHOD{clprnt.end}")
 
    mc_dict = read_MasterControl(input_mcfile)


    ## COLLECT NECESSARY DATA FOR INFERRING THE STARTING STA
    # collect the guide tree from the user or the previous stage
    if input_guide_tree == None:
        guide_tree = get_known_BPP_param(mc_dict, "A00")["newick"]
    else:
        guide_tree = input_guide_tree  
    # collect the imap from the user or the previous stage
    if input_imap == None:
        input_imap = Imap_to_List(get_known_BPP_param(mc_dict, "A00")["Imapfile"])
        indpop_dict = Imap_to_IndPop_Dict(input_imap)
    else:
        indpop_dict = Imap_to_IndPop_Dict(input_imap)
    

    ## SET UP THE STARTING STATE
    HMmode = get_HM_parameters(mc_dict)["mode"]
    accepted_pops, halt_pop_number = get_HM_StartingState(guide_tree, HMmode)
    accepted_pops_over_time = []
    if   HMmode == "merge":
        accepted_pops_over_time.append(accepted_leaves(guide_tree, accepted_pops))
    elif HMmode == "split":
        accepted_pops_over_time.append(accepted_pops)

    ## PRINT INTRODUCTORY FEEDBACK TO THE USER ABOUT THE STARTING STATE
    if   HMmode == "merge":
        print("\nTHE HIERARCHICAL METHOD WILL BE APPLIED IN 'MERGE' MODE.\n")
        print("This means that the program will start with populations separated according to the guide tree,")
        print("and find the optimal delimitation pattern by progressively merging populations.")
    elif HMmode == "split":
        print("\nTHE HIERARCHICAL METHOD WILL BE APPLIED IN 'SPLIT' MODE.\n")
        print("This means that the program will start with all populations merged at the root of the guide tree,")
        print("and find the optimal delimitation pattern by progressively splitting populations.")

    print(f"\n>> STARTING STATE:\n")
    print("THE ACCEPTED SPECIES ARE:\n")
    for population in accepted_pops_over_time[0]:
        print(f"\t{population}")
    print()
    print(f"GUIDE TREE NEWICK:\n\n\t{string_limit(guide_tree, 120)}")
    print("\nGUIDE TREE ASCII:")
    print(tree_ASCII(guide_tree))
    print(f"\nSTARTING IMAP:\n")
    pretty(Imap_to_PopInd_Dict(input_imap))

    ## PERFORM CHECK OF SUITABILITY FOR GDI CALCULATIONS
    print("COMPATIBILITY CHECKING:\n")
    check_GuideTree_Imap_MSA_compat(guide_tree, input_imap, get_known_BPP_param(mc_dict, "A00")["seqfile"])

           # HM ITERATIONS #
    ##############################
    #-----------------------------#
    step = 1
    to_iterate = True
    # run the HM until no more merges or splits can be executed
    while to_iterate == True:
        accepted_pops, to_iterate = HMIteration(input_mcfile        = input_mcfile,
                                                input_guide_tree    = guide_tree,
                                                input_indpop_dict   = indpop_dict,
                                                input_accepted_pops = accepted_pops,
                                                halt_pop_number     = halt_pop_number,
                                                step                = step)
        accepted_pops_over_time.append(accepted_leaves(guide_tree, accepted_pops))
        step += 1
    #-----------------------------#
    ###############################

    ## PROVIDE FEEDBACK TO USER
    print(f"{clprnt.BLUE}\n<< HIERARCHICAL METHOD FINISHED >>{clprnt.end}")
    
    # write final output state to an output folder
    target_dir = f'{input_mcfile[0:-4]}_Final_Result'
    create_TargetDir(target_dir_name=target_dir)
    print("Imap file, phylogenetic tree, and associated images corresponding to final results can be found in the above folder.")
    final_folder = f'{input_mcfile[0:-4]}_2_HM_{step-1}'
    shutil.copy(src = os.path.join(final_folder, "OUTPUT_IMAP.txt"),  dst = target_dir)
    shutil.copy(src = os.path.join(final_folder, "OUTPUT_TREE.txt"), dst = target_dir)
    shutil.copy(src = os.path.join(final_folder, "Output_Topology.png"), dst = target_dir)
    shutil.copy(src = os.path.join(final_folder, "Imap_Visualization.png"), dst = target_dir)
    
    # print the accepted species at each step
    print("\nACCEPTED SPECIES AT EACH ITERATION:\n")
    for i, poplist in enumerate(accepted_pops_over_time):
        if   i == 0:
            print("start:")
            for pop in poplist:
                print(f"\t{pop}")
        elif i+1 == len(accepted_pops_over_time):
            print(f"\nfinal result {i}:")
            for pop in poplist:
                    print(f"\t{pop}")
        else:
            print(f"\niteration {i}:")
            for pop in poplist:
                    print(f"\t{pop}")
    
    print("\n-- end of run --\n")
    exit()