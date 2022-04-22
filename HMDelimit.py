from sys import argv

from check_param_module import check_Master_Control

from controlflow_module import find_initial_State
from controlflow_module import controlled_BPP_cfile_check
from controlflow_module import controlled_BPP_parameter_check

from stage_modules import StartingTopolgy
from stage_modules import StartingDelimitation
from stage_modules import HierarchicalMethod

from cmdline_module import cmdline_interpret

def delimit_steps  (
    mc_file,
    p_state,
                    ):

    if   p_state == 3:
        HierarchicalMethod(mc_file)
        
    elif p_state == 2:
        guide_tree, imap = StartingDelimitation(mc_file)
        HierarchicalMethod(mc_file, input_guide_tree = guide_tree, input_imap = imap)
    
    elif p_state == 1:
        tree = StartingTopolgy(mc_file)
        guide_tree, imap = StartingDelimitation(mc_file, starting_tree = tree)
        HierarchicalMethod(mc_file, input_guide_tree = guide_tree, input_imap = imap)

def HMpipeline(mc_file, checkonly):
    # check if any of the parameters in master control file are erroneously specified
    check_Master_Control(mc_file)
    
    # check which trees are provided, and set the pipeline to begin the appropriate stage
    p_state = find_initial_State(mc_file)
    
    # check if the user supplied BPP control files relevant to the specific steps are correct
    controlled_BPP_cfile_check(mc_file, p_state)
    
    # perform a compatibility check for BPP paramters
    controlled_BPP_parameter_check(mc_file, p_state)

    # exit if pipeline is in check only mode
    if checkonly == True: exit()

    #run the appropriate stages of the pipeline
    delimit_steps(mc_file, p_state)



### ---- MAIN ---- ###    
mc_file, checkonly = cmdline_interpret(argv)
HMpipeline(mc_file, checkonly)