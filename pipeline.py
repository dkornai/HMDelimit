from check_param_module import check_Master_Control
from controlflow_module import find_initial_State

from stage_modules import StartingTopolgy
from stage_modules import StartingDelimitation
from stage_modules import HierarchicalMethod

def HMpipeline(input_control_file):
    # check if any of the parameters in master control file are erroneously specified
    check_Master_Control(input_control_file)

    # # check which trees are provided, and set the pipeline to begin the appropriate stage
    # p_state = find_initial_State(input_control_file)
    
    # tree = StartingTopolgy(input_control_file)
    # guide_tree, imap = StartingDelimitation(input_control_file, starting_tree = tree)
    # HierarchicalMethod(input_control_file, guide_tree, imap)

HMpipeline("micro_test.txt")
HMpipeline("badMC.txt")