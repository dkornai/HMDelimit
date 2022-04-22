'''
THIS MODULE CONTAINS THE FUNCTIONS REQUIRED FOR INTERPRETING THE 
COMMAND LINE ARGUMENTS TO THE PROGRAM.
'''

import os
import os.path
import ntpath
import shutil
from pathlib import Path

from check_helper_functions import check_File_exists, check_Imap_filetype, check_MSA_filetype

from custom_types import Master_control_dict, Master_control_file, file_path

from data_dicts import MCF_param_dict

from helper_functions import pretty

## DATA DEPENDENCIES
from data_dicts import clprnt

## SPECIALIZED HELPER FUNCTIONS

# get just the filename part of a full file path
def path_filename   (
        path:               file_path
                    ):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

## COLLECTION OF ARGUMENTS FROM COMMAND LINE

# collect a list of possible arguments from the command line
def collect_cmdline_args (argument_list):
    # cut off the first argument, as that is just the name of the python file
    argument_list = argument_list[1:]

    # verify that arguments were provided
    if len(argument_list) == 0:
        print("[X] ERROR: PIPELINE CALLED WITHOUT ANY SPECIFEID PARAMETERS")
        exit()
    
    # check that no parameters were provided with spaces between them:
    if "=" in argument_list:
        print("[X] ERROR: ARGUMENT LIST CONTAINS FLOATING '=' VALUES")
        print("Please provide all parameters as: example_parameter1=example_value1 example_parameter2=example_value2")
        exit()
    
    # verify that no duplicate arguments exist
    if len(set(argument_list)) < len(argument_list):
        print("[X] ERROR: THE SAME ARGUMENT IS CALLED MULTIPLE TIMES!")
        exit()

    # see if the "--check" argument is passed, which activates checking only mode
    output_dict = {}
    output_dict["checkonly"] = False
    if "--check" in argument_list:
        output_dict["checkonly"] = True
        argument_list.remove("--check")

    # if one argument is present, check that it is the master control file
    if len(argument_list) == 1:
        # if it is the mcf argument
        if "mcf=" in argument_list[0]:
            try:
                mcf_name = argument_list[0].split("=")[1]
                output_dict["mcf"] = False
                return output_dict

            except:
                print("[X] ERROR: MASTER CONTROL FILE NAME EMPTY")
                exit()
        # if it is not
        else:
            print("[X] ERROR: WHEN ONLY A SINGLE ARGUMENT IS PROVIDED, IT MUST BE A MASTER CONTROL FILE\n")
            print("set the master control file argument as: mcf=example_mastercontrol_file.txt")
            exit()

    # collect the arguments that match the possible parameters of the master control file
    
    output_mc_dict = {}
    correct_arguments = []
    for argument in argument_list:
        for paramname in MCF_param_dict:
            if paramname in argument:
                output_mc_dict[paramname] = argument.split("=")[1]
                correct_arguments.append(argument)
                break
    
    # verify that no arguments supplied remained unmatched
    unrecognized = set(argument_list).difference(set(correct_arguments))
    if len(unrecognized) > 0:
        print("[X] ERROR: THE FOLLOWING ARGUMENTS ARE NOT RECOGNIZED:")
        for param in unrecognized:
            print(param)
            if "mcf=" in param:
                print("WARNING: MASTER CONTROL FILE NAME IS NOT A VALID PARAMETER WHEN SUPPLYING MULTIPLE ARGUMENTS")
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
        print("INITIALIZING PIPELINE WITH USER SUPPLIED MASTER CONTROL FILE!")
        filepath = Path(filepath)
        my_abs_path = filepath.resolve(strict=True)
        
        
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

    filename = "AutoMasterControl.txt"

    print(f"\nAUTOMATICALLY GENERATING MASTER CONTROL FILE '{filename}' BASED ON USER INPUT!\n")

    # check that such a file exists already, and fail the pipeline if so
    file_exists = check_File_exists(filename)
    if file_exists == 1:
        print(f"[X] ERROR: FILE WITH NAME '{filename}' ALREADY EXISTS IN THE DIRECTORY '{os.getcwd()}'")
        exit()

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

    #### FIX VERIFY THAT THE TARGET DIRECTORY IS LEGITIMATE    

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

    param = collect_cmdline_args(argument_list)
    
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