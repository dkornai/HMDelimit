from check_param_module import check_Master_Control

from controlflow_module import controlled_check
from controlflow_module import find_initial_State

from helper_functions import set_wd

from stage_modules import StartingTopolgy
from stage_modules import StartingDelimitation
from stage_modules import HierarchicalMethod

def HMpipeline(input_control_file):
    # start by moving to the correct working directory
    mc_file = set_wd(input_control_file)
    # check if any of the parameters in master control file are erroneously specified
    check_Master_Control(mc_file)
    # check which trees are provided, and set the pipeline to begin the appropriate stage
    p_state = find_initial_State(mc_file)
    # perform a compatibility check for BPP paramters
    compatible = controlled_check(mc_file, p_state)
    # if compatible:
    #     if   p_state == 3:
    #         HierarchicalMethod(mc_file)
        
    #     elif p_state == 2:
    #         guide_tree, imap = StartingDelimitation(mc_file)
    #         HierarchicalMethod(mc_file, guide_tree, imap)
        
    #     elif p_state == 1:
    #         tree = StartingTopolgy(mc_file)
    #         guide_tree, imap = StartingDelimitation(mc_file, starting_tree = tree)
    #         HierarchicalMethod(mc_file, guide_tree, imap)

#HMpipeline("Data_Sets/Rotaria_2009/micro_test.txt")
HMpipeline("Data_Sets/Ants_2021/MC.txt")