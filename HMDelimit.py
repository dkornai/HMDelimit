from sys import argv
from io import StringIO

from check_param_module import check_Master_Control

from controlflow_module import find_initial_State
from controlflow_module import controlled_BPP_cfile_check
from controlflow_module import controlled_BPP_parameter_check

from helper_functions import set_wd
from helper_functions import initialize

from stage_modules import StartingTopolgy
from stage_modules import StartingDelimitation
from stage_modules import HierarchicalMethod

def HMpipeline(mc_file):
    # check if any of the parameters in master control file are erroneously specified
    check_Master_Control(mc_file)
    
    # check which trees are provided, and set the pipeline to begin the appropriate stage
    p_state = find_initial_State(mc_file)
    
    # check if the user supplied BPP control files relevant to the specific steps are correct
    controlled_BPP_cfile_check(mc_file, p_state)
    
    # perform a compatibility check for BPP paramters
    compatible = controlled_BPP_parameter_check(mc_file, p_state)
    
    #run the appropriate stages of the pipeline
    if compatible:
        if   p_state == 3:
            HierarchicalMethod(mc_file)
        
        elif p_state == 2:
            guide_tree, imap = StartingDelimitation(mc_file)
            HierarchicalMethod(mc_file, guide_tree, imap)
        
        elif p_state == 1:
            tree = StartingTopolgy(mc_file)
            guide_tree, imap = StartingDelimitation(mc_file, starting_tree = tree)
            HierarchicalMethod(mc_file, guide_tree, imap)


argument_list = argv
if len(argument_list) == 2:
    master_control_file = argument_list[1]
    print(master_control_file)
    HMpipeline(master_control_file)

elif len(argument_list) == 4:
    master_control_file = initialize(argument_list[1], argument_list[2], argument_list[3])
    HMpipeline(master_control_file)


#"Data_Sets/Rotaria_2009/micro_test.txt"
#"Data_Sets/Ants_2021/MC.txt"