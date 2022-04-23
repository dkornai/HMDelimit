'''
THIS MODULE CONTAINS AN EXTENSIVE ARRAY OF FUNCTIONS THAT AIM
TO CHECK FOR MISSPECIFICATIONS OF PARAMETERS, EITHER IN THE 
MASTER CONTROL FILE, OR BPP CONTROL FILES.
'''
## DEPENDENCDIES
# STANDAR LIBRARY DEPENDENCIES
import os
import copy
import warnings
import re
from pathlib import Path
from custom_types import Master_control_file

# EXTERNAL LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

# HELPER FUNCTION DEPENDENCIES
from helper_functions import alignfile_to_MSA
from helper_functions import readLines
from helper_functions import bppcfile_to_dict
from helper_functions import Imap_to_IndPop_Dict
from helper_functions import Imap_to_List
from helper_functions import Imap_to_PopInd_Dict

## DATA DEPENDENCIES
from data_dicts import HM_decision_criteria
from data_dicts import valid_BPP_param_names
from data_dicts import MCF_param_dict
from data_dicts import default_BPP_param


## CHECK IF THE PARAMETER IS EXPLICITLY WRONGLY SPECIFIED
'''
These functions check if the user has specified a parameter
with the wrong type of value, or a nonsensical value
'''
# check if a supplied filename actually points to an existing file in the current directory
def check_File_exists(path):
    if path == "?":
            file_state = 0
    else:
        try:
            filepath = Path(path)
            my_abs_path = filepath.resolve(strict=True)
            if Path(my_abs_path).is_file():
                file_state = 1
            else:
                file_state = -1
        except:
            file_state = -1

    return file_state

# check if a supplied tree is correctly, formatted, and contains no polytomies
def check_Newick(tree):
    if tree == "?":
        tree_state = 0
    else:
        try:
            # check if the tree can be read in at all
            t = Tree(tree)
            tree_state = 1
            ##### FIX FIX FIX check if the tree contains repeated node names, or conflicts that will arise from internal node names
            
            # check if the tree contrins polytomies, and return a special error if so
            t_2 = copy.deepcopy(t)
            t_2.resolve_polytomy(recursive=True)
            rf = t.compare(t_2)["rf"]
            if rf != 0:
                tree_state = -2

        except:
            tree_state = -1
    
    return tree_state

# check if a single numeric parameter is supplied in the correct format and range
def check_Numeric(value, float_or_int, min = -1000000000, max = 100000000):
    numeric_state = -1
    try:
        num = float(value)
        if float_or_int == "i":
            if num.is_integer() == True and min < num < max:
                numeric_state = 1
        elif float_or_int == "f":
            if min < num < max:
                numeric_state = 1
    except:
        if value == "?":
            numeric_state = 0
    
    return numeric_state

# check if a value is from a presupplied list of acceptable results
def check_ValueIsFrom(input_value, valid_list):
    if input_value == "?":
        value_state = 0
    else:
        value_state = -1
        try:
            if any([input_value == valid for valid in valid_list]):
                value_state = 1
        except:
            value_state = -1

    return value_state

## BPP SPECIFIC MISSPECIFICATION CHECKING FUNCTIONS
##### FIX FIX FIX ADD EXTRA SPECIESMODELPRIOR PARAMETER TO CHECK WHEN CHECKING FOR A11
### FIX FIX FIX ADD LOCUSRATE CHECKER!!!
'''
These functions check if the user has misspelled or
misformatted BPP control file parameters
'''
# check if the finetune parameter passed to BPP is correctly specified
def check_Finetune(finetune):
    if finetune == "?":
        ft_state = 0
    else:
        try:
            ft_state = -1
            ft = finetune.split(" ")
            # check if the first finetune is 0 or 1 followed by a ":", and the remainders are integers
            if float(ft[0][0]) == 0 or float(ft[0][0]) == 1:
                if ft[0][1] == ":":
                    if all(isinstance(float(x), float) for x in ft[1:]) and all(0 < float(x) < 10 for x in ft[1:]):
                        ft_state = 1
                
        except:
            ft_state = -1

    return ft_state

# check if the tau prior parameter passed to BPP is correctly specified
def check_Tauprior(tauprior):
    if tauprior == "?":
        tau_state = 0
    else:
        try:
            tau_state = -1
            t = tauprior.split(" ")
            # check if the tau is two floats, and the second is less than 1
            if len(t) == 2 and 0 < float(t[1]) < 1:
                if all(isinstance(float(x), float) for x in t):
                    tau_state = 1
        except:
            tau_state = -1

    return tau_state

# check if the theta prior parameter passed to BPP is correctly specified
def check_Thetaprior(thetaprior):
    if thetaprior == "?":
        theta_state = 0
    else:
        try:
            theta_state = -1
            t = thetaprior.split(" ")
            # check if the theta is two floats followed by "e", the second is less than 1
            if len(t) == 3 and 0 < float(t[1]) < 1 and t[2] == "e":
                if all(isinstance(float(x), float) for x in t[0:1]):
                    theta_state = 1
        except:
            theta_state = -1

    return theta_state

# check if the threads parameter passed to BPP is correctly specified
'''
This quite complex function checks if the number, starting point, and offset of the 
requested threads is actually compatible with the number of threads available
to the machine. This is beacuse if the user requests more cores than available, or tries to
request that BPP pins threads to nonexistent cores, BPP will crash. 
'''
def check_Threads(threads):
    if threads == "?":
            threads_state = 0
    else:
        # find the available number of CPU cores
        n_cpu = int(os.cpu_count())
        try:
            threads_state = -1
            th = threads.split(" ")
            
            # check if threads is 3 integers, and all values are at least 1 (there is no such thing as 0 threads)
            if all(isinstance(int(x), int) for x in th) and all (int(num) > 0 for num in th):
                th = [int(x) for x in th]
                # if only one integer is given, that is interpreted as the number of threads requested
                if len(th) == 1:
                    # check that the number of threads requested <= threads in the CPU
                    if th[0] <= n_cpu:
                        threads_state = 1
                    else:
                        threads_state = -2
                elif len(th) == 2:
                    # check that the requested offset and the number of threads still fits the CPU
                    if th[0] + (th[1]-1) <= n_cpu:
                        threads_state = 1
                    else:
                        threads_state = -2
                elif len(th) == 3:
                    # check that all the requested threads still fit the CPU
                    if (th[1]-1) + (th[2]*th[0]) <= n_cpu:
                        threads_state = 1
                    else:
                        threads_state = -2
        except:
            threads_state = -1

    return threads_state

# check that the number of threads requested <= the number of loci in the MSA
def check_Threads_MSA_compat(input_threads, alignmentfile):
    n_threads = int(input_threads.split()[0])
    true_nloci = len(alignfile_to_MSA(alignmentfile))

    if n_threads <= true_nloci:
        threads_state = 1
    else:
        threads_state = -3

    return threads_state

# check that the number of threads requested <= the number of loci in the MSA
def check_Threads_nloci_compat(input_threads, input_nloci):
    n_threads = int(input_threads.split()[0])

    if n_threads <= int(input_nloci):
        threads_state = 1
    else:
        threads_state = -4

    return threads_state

# check that the number of loci to check is less than or equal to the loci in the MSA
def check_nloci_MSA_compat(input_nloci, alignmentfile):
    true_nloci = len(alignfile_to_MSA(alignmentfile))
    user_nloci = int(input_nloci)

    if user_nloci <= true_nloci:
        nloci_state = 1
    else:
        nloci_state = -2

    return nloci_state

# check if the print paramters are correctly formatted
def check_Print(printparam):
    if printparam == "?":
        printstate = 0
    else:
        try:
            printsplit = printparam.split(" ")
            if all(char == "0" or char == "1" for char in printsplit):
                printstate = 1
            else:
                printstate = -1
        except:
            printstate = -1

    return printstate

# check that the target output files are single files and not in other folders.
def check_Outfilename(filename):
    if filename == "?":
        file_state = 0
    else:
        try:
            file_state = -1
            if filename.split()[0] == filename and re.match(r"\b[a-zA-Z0-9_]+\.txt$", filename):
                file_state = 1
        except:
            file_state = -1
    
    return file_state

# check that the species&tree line (which also corresponds to my popsizes line) is formatted according to expectations
def check_SandT_popsizes(s_and_t, popsizes):
    if s_and_t == "?" and popsizes == "?":
        st_state = 0
    else:
        st_state = -1
        try:
            st_split = s_and_t.split()
            num = int(st_split[0]) # this should be the number of species
            
            popsize_split = popsizes.split(None); popsize_split = [float(size) for size in popsize_split]
            # if the first value is an integer, and the line contains the required amount of species
            if len(st_split) == num+1:
                # check for repeating population names
                if len(st_split[1:]) == len(set(st_split[1:])):
                    # if the popsizes row also has the correct number of integer values
                    if all(size.is_integer() for size in popsize_split) and len(popsize_split) == num:
                        st_state = 1
        except:
            pass

    return st_state

# check that the parameters given for the speciesdelimitaion line are correctly formatted, if this paramter is supposed to start with 1
def check_speciesdelimitation(sd):
    if sd == "?":
        sd_state = 0
    else:
        try:
            sd_state = -1 
            sdplit = sd.split()
            # if the line is only "1", it is accepted
            if len(sdplit) == 1:
                if sdplit[0] == "1":
                    sd_state = 1
            # if the line starts with "1 0..." this corresponds to eq3 and eq4 in Yang & Rannala (2010)
            elif sdplit[0] == "1" and sdplit[1] == "0" and len(sdplit) == 3:
                if 0 <float(sdplit[2]) < 10:
                    sd_state = 1
            # if the line starts with "1 1..." this corresponds to eq6 and eq7 in Yang & Rannala (2010)
            elif sdplit[0] == "1" and sdplit[1] == "1" and len(sdplit) == 4:
                if 0 <float(sdplit[2]) < 10 and 0 <float(sdplit[3]) < 10:
                    sd_state = 1
        except:
            sd_state = -1 
    
    return sd_state

# check that the speciesdelimitation and speciestree flags match the mode of BPP the control file is intended for
def check_BPP_mode(sd, st, BPP_mode):
    if sd == "?":
        sd_state = 0
    else:
        try:
            sd_state = -1
            if BPP_mode == "A00":
                if   sd == "0":
                    sd_state = 1
                elif check_speciesdelimitation(sd):
                    sd_state = -2
            elif BPP_mode == "A01":
                if   sd == "0":
                    sd_state = 1
                elif check_speciesdelimitation(sd):
                    sd_state = -3
            elif BPP_mode == "A11":
                if   check_speciesdelimitation(sd):
                    sd_state = 1
                elif sd == "0":
                    sd_state = -4
        except:
            sd_state = -1
    
    if st == "?":
        st_state = 0
    else:
        try:
            st_state = -1
            if BPP_mode == "A00":
                if   st == "0":
                    st_state = 1
                elif st == "1":
                    st_state = -2
            elif BPP_mode == "A01":
                if   st == "1":
                    sd_state = 1
                elif st == "0":
                    sd_state = -3
            elif BPP_mode == "A11":
                if   st == "1":
                    sd_state = 1
                elif sd == "0":
                    sd_state = -4
        except:
            st_state = -1
        
    return sd_state, st_state

## HM PARAMETER SPECIFIC CHECKERS
'''
These functions check for simple misspecifications or mispellings of 
HM specific parameters
'''

# check what the merge decisions are based on, and if all parameters are provided to make the requested decisions
def check_GDI_params(gdiparams, param_checked):
    # checks against the list of possible decision criteria (located in "data_dicts")
    possible_answers = HM_decision_criteria
    
    if gdiparams == "?":
        gdiparams_state = 0
    elif gdiparams in possible_answers:
        # if age is one of the decision criteria, but no mutation rate is specified, 
        # this causes a conflict, as node ages in generations cannot be calculated without 
        # knowing the mutation rate
        if gdiparams in ["age","age_&_one_gdi"] and param_checked['mutationrate'] != 1:
            gdiparams_state = -2
        else:
            gdiparams_state = 1
    else:
        gdiparams_state = -1

    return gdiparams_state


## CHECK THE TYPE AND INTEGRITY OF SPECIFIC SUPPLIED FILES
# check that the Master Control file is valid, and only contains valid parameters
'''
This function aims to ensure that the master control file only contains calls
to valid parameters of the pipeline. In other cases, the pipeline may fail in
unexpected ways, which should be avoided as much as possible.
'''
def check_Master_control_filetype(input_mc_file):
    ## CHECK IF A FILE IS SUPPLIED, OR EXISTS IN THE FIRST PLACE
    mcf_status = check_File_exists(input_mc_file)
    if mcf_status == 0:
        print("[X] ERROR: NO MASTER CONTROL FILE SPECIFIED!")
        print("\t\t\n-- EXITING PROGRAM --")
        exit()
    elif mcf_status == -1:
        print(f"[X] ERROR: NO FILE OF ANY TYPE AT THE REQUESTED LOCATION: {input_mc_file})")
        print("\t\t\n-- EXITING PROGRAM --")
        exit()
    
    ## CHECK IF ALL LINES OF THE MASTER CONTROL FILE CORRESPOND TO CORRECTLY NAMED KNOWN PARAMETERS
    elif mcf_status:
        # load file, and read into numbered rows
        lines = readLines(input_mc_file)
        lines_num = {i:line for i, line in enumerate(lines)}
            # strip comments, and empty lines
        lines_num = {key:lines_num[key].split("#")[0] for key in lines_num}
        lines_num = {key:lines_num[key].split("*")[0] for key in lines_num}
        lines_num = {key:lines_num[key] for key in lines_num if re.search("\S+", lines_num[key])}

        matched_lines = {}
        unmatched_lines = copy.deepcopy(lines_num)
        used_params = []
        # check if the remaining lines contain all valid parameter names, and valid parameters are mentioned only once
        for lineindex in lines_num:
            line = lines_num[lineindex]
            for param in MCF_param_dict:
                if param not in used_params:
                    param_keyphrase = MCF_param_dict[param] # get the keyphrase from the param dict
                    if param_keyphrase in line:
                        matched_lines[lineindex] = line
                        del unmatched_lines[lineindex]
                        used_params.append(param)
                        break

        
        # if all lines were used to specify valid MCF paramters, the program can continue
        if set(matched_lines) == set(lines_num):
            print()
        
        # if no matching parameters were found, the file is not a master control file
        elif len(matched_lines) == 0:
            print(f"[X] ERROR: THE FILE '{input_mc_file}' CONTAINS NO MASTER CONTROL PARAMETERS")
            exit()
        
        # if some of the lines are not valid parameters, the file cannot be correctly interpreted
        else:
            # print the lines that were incorrectly found
            print(f"[X] ERROR: THE MASTER CONTROL FILE '{input_mc_file}' CONTAINS UNRECOGNIZED OR DUPLICATE PARAMETERS")
            print("\tTHE FOLLOWING LINES ARE UNRECOGNIZED OR DUPLICATES:\n")
            for nomatch in unmatched_lines:
                print(f"\tLine {nomatch+1}: {unmatched_lines[nomatch]}")
            exit()


# check if file is a valid BPP control file
def check_BPP_ctl_filetype(bpp_ctl_file):
    cfile_state = check_File_exists(bpp_ctl_file)
    if cfile_state == 1:
        try:
            bpp_cdict = bppcfile_to_dict(bpp_ctl_file)
            cfile_state = 1
        except:
            cfile_state = -2 # the file was not able to be interpreted as a BPP command file
    
    return cfile_state


# check that the BPP control file only contains parameters that are valid for BPP
'''
This function aims to check, in advance, if the user supplied BPP control files are valid.
The function checks to see if all parameters of the control file are in the known valid list,
and will fail the check if not. Additionally, the module also warns the user if they
supply a parameter that is ignored by the pipeline. This is needed, because in order to ensure
no crashes, the pipeline only interprets a limited list of ~20 parameters, for which it
performs specific checks.
'''
def check_BPP_ctl_validity(bpp_ctl_file):
    bpp_cdict = bppcfile_to_dict(bpp_ctl_file)

    matched_lines = {}
    ignored_lines = {}
    unmatched_lines = {}
    # check if the remaining lines contain all valid parameter names, and valid parameters are mentioned only once
    for param in bpp_cdict:
        # ignore these parameters, as they are only internal to the pipeline
        if param == "newick" or param == "popsizes":
            continue
        # check if parameter is in the list of valid ones
        if param in valid_BPP_param_names:
            # check if the parameter is in the list of tracked parameters
            if param in default_BPP_param:
                matched_lines[param] = bpp_cdict[param]
            # collect parameters that are valid but not passed
            else:
                ignored_lines[param] = bpp_cdict[param]
        else:
            unmatched_lines[param] = bpp_cdict[param]

    # warn the user about ignored lines
    if len(ignored_lines) > 0:
        print("\n\tWARNING: The following parameters are present in the BPP control file, but will be ignored by the pipeline:\n")
        for param in ignored_lines:
            print(f"\t{param}")
    
    # raise error if unmatched parameters are found
    if len(unmatched_lines) > 0:
        print("\n\t[X] ERROR: THE FOLLOWING PARAMETERS ARE UNKNOWN TO BPP:\n")
        for param in unmatched_lines:
            print(f"\t{param}")

        compat = False
    
    # pass completely if no errors are found
    elif len(unmatched_lines) == 0:
        print(f"\n\t[*] All parameter names in the BPP control file '{bpp_ctl_file}' are valid\n")

        compat = True

    return compat

# check if an alignment file can be loaded in as a valid MSA object
def check_MSA_filetype(alignmentfile):
    align_state = check_File_exists(alignmentfile)
    if align_state == 1:
        try:
            # try to load the alignment file to the internal MSA object
            align = alignfile_to_MSA(alignmentfile)
            align_state = 1
        except:
            align_state = -2 # the file could not be interpreted as a phylip MSA

    return align_state

# check if the file supposted to be an imap is actually an Imap
def check_Imap_filetype(imapfile):
    imap_state = check_File_exists(imapfile)
    if imap_state == 1:
        try:
            # try to load the imap in all three modes
            imap = Imap_to_List(imapfile)
            imap = Imap_to_PopInd_Dict(imapfile)
            imap = Imap_to_IndPop_Dict(imapfile)
            imap_state = 1
        except:
            imap_state = -2 # lower errors than -1 correspond to specific problems

    return imap_state

# check that the target output folders of the pipeline do not exist
def check_folders_do_not_exist  (
        input_mcfile:                   Master_control_file
                                ):

    dir_list = os.listdir()

    outfolder = []
    outfolder.append(f'{input_mcfile[0:-4]}_1_StartDelim')
    outfolder.append(f'{input_mcfile[0:-4]}_0_StartPhylo')
    outfolder.append(f'{input_mcfile[0:-4]}_Final_Result')
    for i in range(1, 51):
        outfolder.append(f'{input_mcfile[0:-4]}_2_HM_{i}')

    if any(folder in dir_list for folder in outfolder):
        print(f"[X] ERROR: OUTPUT FOLDERS FOR {input_mcfile} ALREADY PRESENT IN THE DIRECTORY!")
        print("THE FOLLOWING FOLDERS ARE CAUSING THE ERROR:")
        for element in dir_list:
            if element in outfolder:
                print(element)
        print("\nConsider changing the working directory, master control file name, or moving previous results to a new folder")
        exit()