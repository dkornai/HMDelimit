'''
THIS MODULE CONTAINS THE FUNCTIONS REQUIRED FOR INTERPRETING THE 
COMMAND LINE ARGUMENTS TO THE PROGRAM.
'''
## DEPENDENCIES
# STANDARD LIBRARY
import os
import os.path
import ntpath
import shutil
from pathlib import Path
import difflib

# CHECK HELPERS
from check_helper_functions import check_File_exists
from check_helper_functions import check_Imap_filetype
from check_helper_functions import check_MSA_filetype

# HELPER FUNCTIONS
from helper_functions import pretty
from helper_functions import stripall

## TYPE HINTS
from custom_types import Master_control_dict
from custom_types import Master_control_file
from custom_types import file_path

## DATA DEPENDENCIES
from data_dicts import clprnt
from data_dicts import MCF_param_dict
from data_dicts import command_line_params


## SPECIALIZED HELPER FUNCTIONS

# get just the filename part of a full file path
def path_filename   (
        path:               file_path
                    ):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

## GETTING PIPELINE ARGUMENTS FROM THE COMMAND LINE
# read in the raw command line parameters, and ensure they are formatted correctly.
def collect_cmdline_args    (
        argument_list:              list[str]
                            ) ->    list[str]:
    
    # cut off the first argument, as that is just the name of the python file
    argument_list = argument_list[1:]

    # verify that arguments were provided
    if len(argument_list) == 0:
        print("[X] ERROR: PIPELINE CALLED WITHOUT ANY SPECIFEID PARAMETERS")
        exit()
    
    # check that parameters are separated by commas
    argument_string = " ".join(argument_list)
    if "," not in argument_string:
        print("[X] ERROR: ARGUMENTS MUST BE SEPARATED BY ','")
        print("Please provide all parameters as: 'example_parameter1 = example_value1, example_parameter2 = example_value2'")
        print("even when providing only the 'mcf' parameter, provide it as: 'mcf=example_mcf,'")
        exit()
    
    # paste the arguments together,cut up at ',', remove extra whitespaces, and remove any components that are now empty
    argument_list = argument_string.split(',')
    argument_list = [stripall(argument) for argument in argument_list]
    argument_list = [argument for argument in argument_list if len(argument) > 0]

    # verify that no duplicate arguments exist
    if len(set(argument_list)) < len(argument_list):
        print("[X] ERROR: THE SAME ARGUMENT IS CALLED MULTIPLE TIMES!")
        exit()

    return argument_list

# interpret the incoming command line arguments
def interpret_cmdline_args  (
        argument_list:              list
                            ) ->    dict:
    
    # filter and chcek for misspecified arguments
    argument_list = collect_cmdline_args(argument_list)

    # see if the "--check" argument is passed, which activates checking only mode
    output_dict = {}
    output_dict["checkonly"] = False
    if "--check" in argument_list:
        output_dict["checkonly"] = True
        argument_list.remove("--check")

    # if the master control file argument is present, check that it is the only one
    
    # if it is the mcf argument
    if "mcf" in argument_list[0] and len(argument_list) ==1:
        try:
            mcf_name = argument_list[0].split("=")[1]
            mcf_name = stripall(mcf_name)
            output_dict["mcf"] = stripall(mcf_name)
            return output_dict

        except:
            print("[X] ERROR: MASTER CONTROL FILE NAME EMPTY")
            exit()
    elif "mcf" in argument_list[0]:
        print("[X] ERROR: WHEN MASTER CONTROL FILE ARGUMENT IS PROVIDED, NO OTHER PARAMETERS CAN BE SET FROM THE COMMAND LINE")
        exit()
    # if it is not
    elif len(argument_list) == 1:
        print("[X] ERROR: WHEN ONLY A SINGLE ARGUMENT IS PROVIDED, IT MUST BE A MASTER CONTROL FILE\n")
        print("set the master control file argument as: mcf=example_mastercontrol_file.txt")
        exit()

    # collect the arguments that match the possible parameters of the master control file
        # parameters are the standard MCF parameters, + an additional working directory parameter
    pipeline_params = command_line_params
    output_mc_dict = {}
    correct_arguments = []
    for argument in argument_list:
        for paramname in pipeline_params:
            if paramname in argument:
                output_mc_dict[paramname] = stripall(argument.split("=")[1])
                correct_arguments.append(argument)
                break
    
    # verify that no arguments supplied remained unmatched
    unrecognized = set(argument_list).difference(set(correct_arguments))
    if len(unrecognized) > 0:
        print("[X] ERROR: THE FOLLOWING ARGUMENTS ARE NOT RECOGNIZED:")
        for param in unrecognized:
            # feedback to user if a close match is found
            closest_match = difflib.get_close_matches(param, pipeline_params, 1, 0.5)
            match_text = ''
            if len(closest_match) > 0:
                match_text += f"\t -- did you mean '{str(closest_match)[2:-2]}'?"
            print(f"\t{param}{match_text}")
            # feedback to user if they tried to supply a master control file among multiple arguments
            if "mcf=" in param:
                print("WARNING: MASTER CONTROL FILE NAME IS NOT A VALID PARAMETER WHEN SUPPLYING MULTIPLE ARGUMENTS")
            if "tree" in param:
                print("WARNING: TREE PARAMETERS CAN NOT BE PROVIDED FROM THE COMMAND LINE, ONLY THE MASTER CONTROL FILE")
        exit()
    
    # verify that no matched arguments contain empty values
    zerolength = [param for param in output_mc_dict if len(output_mc_dict[param]) == 0]
    if len(zerolength) > 0:
        print("[X] ERROR: THE FOLLOWING PARAMETERS ARE PROVIDED WITH EMPTY VALUES:")
        for param in zerolength:
            print(param)
        exit()
            

    # check that the three minimum arguments were present in the call
    if any(param not in output_mc_dict for param in ["Imapfile", "seqfile", "working_dir"]):
        print("[X] ERROR: REQUIRED PARAMETERS ARE MISSING!\n")
        print("The pipeline is missing the following parameter:")
        print(str([param for param in ["Imapfile", "seqfile", "working_dir"] if param not in output_mc_dict])[1:-1])
        exit()

    else:
        
        output_dict = output_dict|output_mc_dict #merge master control dict and check only status
        return output_dict


## SET THE CORRECT WORKING DIRECTORY IF THE PIPELINE WAS PROVIDED A MASTER CONTROL FILE
def move_to_mc_folder   (
        mc_name:                dict
                        ) ->    Master_control_file:

    filepath = mc_name["mcf"]

    try:
        print("INITIALIZING PIPELINE WITH USER SUPPLIED MASTER CONTROL FILE")
        filepath = Path(filepath)
        my_abs_path = filepath.resolve(strict=True)
        
        # move the program to working directory
        working_dir = os.path.dirname(my_abs_path)
        os.chdir(working_dir)
        print(f"\nPipeline working directory is: {working_dir}\n")

        return path_filename(my_abs_path)

    except:
        return filepath
    

## AUTOMATIC MASTER CONTROL FILE CREATION

# create a master control file with the arguments provided in the command line
def create_auto_MC  (
        mc_dict:            Master_control_dict,
                    ) ->    Master_control_file:

    filename = "AutoMC.txt"

    print(f"\nAUTOMATICALLY GENERATING MASTER CONTROL FILE '{filename}' BASED ON USER INPUT")

    # check that such a file exists already, and fail the pipeline if so
    file_exists = check_File_exists(filename)
    if file_exists == 1:
        print(f"[X] ERROR: FILE WITH NAME '{filename}' ALREADY EXISTS IN THE DIRECTORY '{os.getcwd()}'")
        exit()

    # strip the working directory parameter, as that is not passed
    mc_dict.pop("working_dir")

    # generate the control file
    f = open(filename, "x")
    parameter_lines = [f"{MCF_param_dict[param]} = {mc_dict[param]}\n" for param in mc_dict]
    f.writelines(parameter_lines)

    return filename

# automatically set up the master control file if called from the command line
def initialize_with_autoMC  (
        mc_dict:                    Master_control_dict,
                            )  ->   Master_control_file:

    # print feedback to the user
    print("INITIALIZING AUTOMATIC MASTER CONTROL FILE CREATION\n")
    pretty(mc_dict)

    seqfile         = mc_dict["seqfile"]
    Imapfile        = mc_dict["Imapfile"]
    target_dir_name = mc_dict["working_dir"]

    # verify that the seqfile exists and is the correct type
    file_verif = {}
    if os.path.isfile(seqfile):
        if check_MSA_filetype(seqfile):
            file_verif["seq"] = True
        else:
            print(f"[X] ERROR: '{seqfile}' IS NOT A VALID PHYLIP MSA")
            file_verif["seq"] = False
    else:
        print(f"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION: {seqfile}")
        file_verif["seq"] = False
    
    # verify that the Imapfile exists and is the correct type
    if os.path.isfile(Imapfile):
        if check_Imap_filetype(Imapfile):
            file_verif["Imap"] = True
        else:
            print(f"[X] ERROR: '{Imapfile}' IS NOT A VALID IMAP FILE")
            file_verif["seq"] = False
    else:
        print(f"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION: {Imapfile}")
        file_verif["Imap"] = False

    ### TO DO: VERIFY THAT THE TARGET DIRECTORY IS LEGITIMATE    

    # exit if any of the file checks failed
    if file_verif["Imap"] == False or file_verif["seq"] == False:
        exit()

    # start my creating the target directory if it does not exist
    if not os.path.exists(target_dir_name):
        os.makedirs(target_dir_name)
        print(f"Working directory created: {target_dir_name}\n")
    else:
        print(f"Existing folder is used as working directory: {target_dir_name}\n")

    # collect just the name part of the files, and see if they exist in the target directory
    file_present ={}
    seqfile_name = path_filename(seqfile)
    mc_dict["seqfile"] = seqfile_name #change this value in the MCF, as the file will be in the working directory
    if Path(os.path.join(target_dir_name, seqfile_name)).is_file():
        file_present["seq"] = True
        print("Specified seqfile already in working directory")
    else:
        file_present["seq"] = False
    
    Imapfile_name = path_filename(Imapfile)
    mc_dict["Imapfile"] = Imapfile_name #change this value in the MCF, as the file will be in the working directory
    if Path(os.path.join(target_dir_name, Imapfile_name)).is_file():
        file_present["Imap"] = True
        print("Specified Imapfile already in working directory")
    else:
        file_present["Imap"] = False

    # if the files are not already in the directory, copy them in
    if file_present["seq"] == False:
        shutil.copy(src = seqfile,  dst = target_dir_name)
        print("seqfile copied to working directory")
    if file_present["Imap"] == False:
        shutil.copy(src = Imapfile, dst = target_dir_name)
        print("Imapfile copied to working directory")


    # move to the target directory
    os.chdir(target_dir_name)

    # create a small master control file: 
    filename = create_auto_MC(mc_dict)

    return filename

## FINAL WRAPPER FUNCTION
def cmdline_interpret   (
    argument_list
                        ):
    print(f"{clprnt.BLUE}<< STARTING HMDELIMIT PIPELINE >>{clprnt.end}\n")

    param = interpret_cmdline_args(argument_list)
    
    # separate out the check only parameter
    checkonly = param["checkonly"]
    param.pop("checkonly")
    if checkonly == True:
        print("\tCHECK ONLY MODE ACTIVATED!\n")

    # this is the case where the pipeline was initalized and pointed to a master control file
    if "mcf" in param:
        mc_file = move_to_mc_folder(param)
    
    # this is the case where the pipeline was initalized with command line specified master control file parameters
    else:
        mc_file = initialize_with_autoMC(param)

    return mc_file, checkonly