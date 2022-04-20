## DEPENDENCDIES
# STANDAR LIBRARY DEPENDENCIES
import os
import copy
import warnings
import re

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

# ALIGNMENT AND IMAP SPECIFIC DEPENDENCIES
from align_imap_module import autoPopParam, count_Seq_Per_Pop

## DATA DEPENDENCIES
from data_dicts import HM_decision_criteria
from data_dicts import valid_BPP_param_names
from data_dicts import MCF_param_dict


## CHECK IF THE PARAMETER IS EXPLICITLY WRONGLY SPECIFIED
'''
These functions check if the user has specified a parameter
with the wrong type of value, or a nonsensical value
'''
# check if a supplied filename actually points to an existing destination
def check_File(filepath):
    file_state = -1
    if filepath == "?":
            file_state = 0
    else:
        try:
            file_exists = os.path.exists(filepath)
            if file_exists == True:
                file_state = 1
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
            # check if the tree contains repeated node names, or conflicts that will arise from internal node names
            ##### FIX FIX FIX 
            
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

## CHECK THE TYPE AND INTEGRITY OF SPECIFIC SUPPLIED FILES
# check that the Master Control file is valid, and only contains valid parameters
'''
This function aims to ensure that the master control file only contains calls
to valid parameters of the pipeline. In other cases, the pipeline may fail in
unexpected ways, which should be avoided as much as possible.
'''
def check_Master_control_filetype(input_mc_file):
    ## CHECK IF A FILE IS SUPPLIED, OR EXISTS IN THE FIRST PLACE
    mcf_status = check_File(input_mc_file)
    if mcf_status == 0:
        print("ERROR: NO MASTER CONTROL FILE SPECIFIED!")
        print("\t\t\n-- EXITING PROGRAM --")
        exit()
    elif mcf_status == -1:
        print(f"ERROR: NO FILE OF ANY TYPE AT THE REQUESTED LOCATION: {input_mc_file})")
        print("\t\t\n-- EXITING PROGRAM --")
        exit()
    
    ## CHECK IF ALL LINES OF THE MASTER CONTROL FILE CORRESPOND TO CORRECTLY NAMED KNOWN PARAMETERS
    elif mcf_status == 1:
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
            print("[*] MASTER CONTROL FILE CALLS ONLY KNOWN PARAMETERS\n")
        
        # if no matching parameters were found, the file is not a master control file
        elif len(matched_lines) == 0:
            print(f"[X] ERROR: THE FILE '{input_mc_file}' CONTAINS NO MASTER CONTROL PARAMETERS")
            print("\t\t\n-- EXITING PROGRAM --")
            exit()
        
        # if some of the lines are not valid parameters, the file cannot be correctly interpreted
        else:
            # print the lines that were incorrectly found
            print(f"[X] ERROR: THE MASTER CONTROL FILE '{input_mc_file}' CONTAINS UNRECOGNIZED OR DUPLICATE PARAMETERS")
            print("\tTHE FOLLOWING LINES ARE UNRECOGNIZED OR DUPLICATES:\n")
            for nomatch in unmatched_lines:
                print(f"\tLine {nomatch+1}: {unmatched_lines[nomatch]}")
            print("\t\t\n-- EXITING PROGRAM --")
            exit()

# check if an alignment file can be loaded in as a valid MSA object
def check_MSA_filetype(alignmentfile):
    align_state = check_File(alignmentfile)
    if align_state == 1:
        try:
            # try to load the alignment file to the internal MSA object
            align = alignfile_to_MSA(alignmentfile)
                # check that the alignment has one loci, and at least two sequences at each loci
            align_state = 1
        except:
            align_state = -2 # lower errors than -1 correspond to specific problems

    return align_state

# check if the file supposted to be an imap is actually an Imap
def check_Imap_filetype(imapfile):
    imap_state = check_File(imapfile)
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

# check if file is a valid BPP control file
def check_BPP_ctl_filetype(bpp_ctl_file):
    cfile_state = check_File(bpp_ctl_file)
    if cfile_state == 1:
        try:
            # attempt to load file as a command dict
            dict = bppcfile_to_dict(bpp_ctl_file)
            
            # check if the command dict only includes parameters that are in the list of known valid ones
            param = set(dict)
            if param.issubset(set(valid_BPP_param_names)):
                cfile_state = 1
            else:
                cfile_state = -3
        except:
            cfile_state = -2 # lower errors than -1 correspond to specific problems
    
    return cfile_state


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
def check_Threads_nloci_compat(input_threads, input_nloci):
    n_threads = int(input_threads.split()[0])
    n_loci = int(input_nloci)

    if n_threads <= n_loci:
        threads_state = 1
    else:
        threads_state = -3

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


## SPECIALIZED DEEP CONFLICT CHECKING FUNCTIONS
'''
These functions exist to check for conflicts between files that are of the correct type.
For example, an Imap file and Alignment file might not be compatible, because individual IDs
exist in the alignment that do not exist in the Imap. 
'''
# check if an Imap file and a sequence alignment are mutually compatible

    # this implies that all individual IDs in the alignment are mapped to a population in the IMAP

def assert_Imap_Seq_compat(imapfile, alignmentfile):
    # get the list of individual IDs in the Imap
    names_imap = set(Imap_to_List(imapfile)[0])
    
    # get the list of individual IDs in the alignment
    alignment = alignfile_to_MSA(alignmentfile)
    names_align = []
    for locus in alignment:
        for seq in locus:
            name = seq.id
            name = name.split("^")[-1]
            names_align.append(name)
    names_align = set(names_align)
    
    # check if the two sets of names are identical
    if names_imap == names_align:
        print("\t[*] The supplied Imap and Sequence alignment are compatible!\n")
        comp_status = 1
    
    # if not, provide detailed feedback about the missing populations
    else:
        print("\t[X] ERROR: THE SUPPLIED IMAP AND SEQUENCE ALIGNMENT ARE INCOMPATIBLE!\n")
        
        not_found_in_alignment = names_imap.difference(names_align)
        if len(not_found_in_alignment) > 0:
            print("\tTHE FOLLOWING INDIVIDUALS FROM THE IMAP ARE NOT FOUND IN THE ALIGNMENT:\n")
            print(f"\t{not_found_in_alignment}\n")
        
        not_found_in_imap = names_align.difference(names_imap)
        if len(not_found_in_imap) > 0:
            print("\tTHE FOLLOWING INDIVIDUALS FROM THE ALIGNMENT ARE NOT FOUND IN THE IMAP:\n")
            print(f"\t{not_found_in_imap}\n")
        comp_status = -1
    
    return comp_status

# wrap the checker into a conditional form that will not run if some of the required parameters are missing
def check_Imap_Seq_compat(param_check, imapfile, alignmentfile):
    # only attempt checking if these files exist in the first place
    if param_check["Imapfile"] != 0 and param_check["seqfile"] != 0:
        print("\tChecking compatibility of user specified Imap file and MSA:\n")
        print(f"\t1) Imap file = {imapfile}")
        print(f"\t2) MSA file  = {alignmentfile}\n")
        
        # check can only proceeed if both files are of the correct type to begin with
        if param_check["Imapfile"] and param_check["seqfile"]:
            comp_status = assert_Imap_Seq_compat(imapfile, alignmentfile)

        else: # if the files exist but have errors, they cannot be checked for compatibility
            print("\t[X] ERROR: MSA-IMAP COMPATIBILTIY CANNOT BE ASSESED DUE TO UPSTREAM MISSPECIFICATIONS") 
            comp_status = -1

    else: # if both or one of the files are missing, dont assess the validity
        comp_status = 0

    return comp_status


# check if a Newick tree and an Imap file are mutually compatible
'''
This implies that the all of the tree leaf names in the newick string are also
found in the Imap file, and vice versa.
'''
def assert_Imap_Tree_compat(imapfile, tree):
    # get the list of population names in the Imap
    pops_imap = set(Imap_to_List(imapfile)[1])
    
    # get the list of populations mentioned in the tree
    pops_tree = []
    t = Tree(tree)
    for node in t.traverse("levelorder"):
        if node.name != "":
            pops_tree.append(node.name)
    pops_tree = set(pops_tree)
    
    # check if the two sets of names are identical
    if pops_imap == pops_tree:
        print("\t[*] The supplied Imap and Newick Tree are compatible!\n")
        comp_status = 1
    
    # if not, provide detailed feedback about the missing populations
    else:
        print("\t[X] ERROR: THE SUPPLIED IMAP AND NEWICK TREE ARE INCOMPATIBLE!\n")
        
        not_found_in_tree = pops_imap.difference(pops_tree)
        if len(not_found_in_tree) > 0:
            print("\tTHE FOLLOWING POPULATIONS FROM THE IMAP ARE NOT FOUND ON THE TREE:\n")
            print(f"\t{not_found_in_tree}\n")
        
        not_found_in_imap = pops_tree.difference(pops_imap)
        if len(not_found_in_imap) > 0:
            print("\tTHE FOLLOWING POPULATIONS FROM THE TREE ARE NOT FOUND IN THE IMAP:\n")
            print(f"\t{not_found_in_imap}\n")
        comp_status = -1   

    return comp_status

# wrap the checker into a conditional form that will not run if some of the required parameters are missing
def check_Imap_Tree_compat(param_check, imapfile, tree):
    # only attempt checking if the paramteters are supplied
    if param_check["newick"] != 0 and param_check["Imapfile"] != 0:
        print("\tChecking compatibility of user specified Imap file and newick tree:\n")
        print(f"\t1) Imap file   = {imapfile}")
        print(f"\t2) Newick tree = {tree}\n")

        # check can only proceed if both parameters are of the correct type to begin with
        if param_check["newick"] and param_check["Imapfile"]:
            comp_status = assert_Imap_Tree_compat(imapfile, tree)
        
        else: # if the files exist but have errors, they cannot be checked for compatibility
            print("\t[X] ERROR: TREE-IMAP COMPATIBILTIY CANNOT BE ASSESED DUE TO UPSTREAM MISSPECIFICATIONS") 
            comp_status = -1
    else: # if both or one of the files are missing, dont assess the validity
        comp_status = 0
    
    return comp_status


# check if the species&tree row implies the same populations structure that the Imap and the seq alignment do.
'''
This implies that the number of sequences that are associated with a population in the data (which inferrable 
from the imap and the MSA, which together tell us which seq IDs in the MSA belong to that population), 
and the species&tree row of the Imap file are identical. If this is not the case, it will cause a conflict
when the pipeline is run (because BPP will search for populations or sequences not in the data), 
so it is crucial that this is checked.
'''
def assert_SandT_Imap_MSA_compat(s_and_t, popsizes, imapfile, alignmentfile):
    # move the user supplied parameters into a dict
    user_parameters = {"species&tree": s_and_t, "popsizes": popsizes }
    
    # generate the same paramters from the data
    generated_parameters = autoPopParam(imapfile, alignmentfile)

    # small function to generate a dict where population names are associated with their sizes
    def deconstruct_popparam(input_dict):
        pop_names = input_dict["species&tree"].split()[1:]
        popsizes = input_dict["popsizes"].split()

        return {pop_names[i]:popsizes[i] for i in range(len(pop_names))}
    
    user_popparam_dict = deconstruct_popparam(user_parameters)
    auto_popparam_dict = deconstruct_popparam(generated_parameters)

    user_pops = set(user_popparam_dict)
    auto_pops = set(auto_popparam_dict)

    # start by checking that the populations referenced in the two files are are identical
    if user_pops != auto_pops:
        print("\t[X] ERROR: THE POPULATIONS REFERENCED IN 'SPECIES&TREE' DO NOT MATCH WITH THE IMAP!")
        not_found_in_cfile = auto_pops.difference(user_pops)
        if len(not_found_in_cfile) > 0:
            print("\tTHE FOLLOWING POPULATIONS FROM THE IMAP ARE NOT FOUND IN THE BPP CFILE:\n")
            print(f"\t{not_found_in_cfile}\n")
        not_found_in_imap = user_pops.difference(auto_pops)
        if len(not_found_in_imap) > 0:
            print("\tTHE FOLLOWING POPULATIONS FROM THE BPP CFILE ARE NOT FOUND IN THE IMAP:\n")
            print(f"\t{not_found_in_imap}\n")
        comp_status = -1
    
    # if they were identical, check if they have identical values for their population sizes
    else:    
        if any(user_popparam_dict[user_pop] != auto_popparam_dict[user_pop] for user_pop in user_pops):
            print("\t[X] ERROR: THE # OF SEQUENCES ASSIGNED TO POPULATIONS IN 'SPECIES&TREE' DO NOT MATCH WITH THE MSA!")
            print("\tTHE FOLLOWING POPULATIONS ARE WRONGLY LABELLED:\n")
            for pop in user_pops:
                if user_popparam_dict[pop] != auto_popparam_dict[pop]:
                    print(f"\t population name: {pop} | # of associated sequences found in the MSA: {auto_popparam_dict[pop]} | # suggested by BPP: {user_popparam_dict[pop]}")
            comp_status = -1
        else:
            print("\t[*] The population names and sizes in 'species&tree' match with those found in the Imap and the MSA!\n")
            comp_status = 1

    return comp_status

# wrap the checker into a conditional form that will not run if some of the required parameters are missing
def check_SandT_Imap_MSA_compat(param_check, s_and_t, popsizes, imapfile, alignmentfile):
    # only attempt the check if all of the parameters are supplied to begin with
    if param_check["species&tree"] != 0 and param_check["Imapfile"] != 0 and param_check["seqfile"] != 0:
        print("\tChecking compatibility of 'species&tree', Imap file, and MSA:\n")
        print(f"\t1) Imap file   = {imapfile}")
        print(f"\t2) MSA file    = {alignmentfile}")
        print(f"\t3) species&tree line supplied in BPP control file:\n")
        print(f"\t\t{s_and_t}")
        print(f"\t\t   {popsizes}\n")
        
        # check can only proceed if all parameters are of the correct type to begin with, and no upstream incompatibilities exist
        if param_check["species&tree"] and param_check["seqfile"] and param_check["Imapfile"] and param_check['imap_seq_compat']:
            comp_status = assert_SandT_Imap_MSA_compat(s_and_t, popsizes, imapfile, alignmentfile)

        else: # if the files exist but have errors, they cannot be checked for compatibility
            print("\t[X] ERROR: SPECIES&TREE-IMAP-MSA COMPATIBILTIY CANNOT BE ASSESED DUE TO UPSTREAM MISSPECIFICATIONS") 
            comp_status = -1
        
    else: # if both or one of the files are missing, dont assess the validity
        comp_status = 0
    
    return comp_status

### FIX FIX INTEGRATE THIS CHECK WHEN STARTING A11 OR A00 calculations
# check if a guide tree and Imap and alignment together are suitable for GDI calculations
'''
This means that leaf node pairs have at least 3 individual samples spread across two populations.
This is required because the GDI relies on the Theta parameters of atleast one of the populations 
that are under consideration for a merge or split. In other cases, the GDI cannot be calculated, 
and the pipeline cannot run, except if it only uses ages.
'''
def check_GuideTree_Imap_MSA_compat(input_newick, imapfile, alignmentfile):
    print("\tChecking if the guide tree, Imap file, and MSA are suitable for GDI calculations:\n")
    print(f"\t1) Imap file   = {imapfile}")
    print(f"\t2) MSA file    = {alignmentfile}")
    print(f"\t3) Guide Tree  = {input_newick}\n")

    popind_dict = Imap_to_PopInd_Dict(imapfile)   
    alignment = alignfile_to_MSA(alignmentfile)

    # count the max numer of sequences per population
    maxcounts = count_Seq_Per_Pop(popind_dict, alignment)

    # collect the population pairs where the number of sequences needs to be checked
    tree = Tree(input_newick)
    pairs_to_test = []
    for node in tree.traverse("levelorder"):
        if len(list(node.iter_descendants("levelorder"))) == 2:
            pair = [leafnode.name for leafnode in node.iter_descendants("levelorder")]
            pairs_to_test.append(pair)
    
    # count the number of combined sequnces in a population pair
    paircounts = {str(pair):(maxcounts[pair[0]] + maxcounts[pair[1]]) for pair in pairs_to_test}
    
    # check for unsuitable pairs
    unsuitable_pairs = [f"\tpair: {pair} \t# of seqs: {paircounts[pair]}" for pair in paircounts if paircounts[pair] < 3]

    # provide detailed feedback if incompatibility is detected
    if len(unsuitable_pairs) > 0:
        print("\t[X] ERROR: THE GUIDE TREE, IMAP, AND SEQUENCE ALIGNMENT ARE UNSUITABLE FOR GDI CALCULATIONS!\n")
        print("\tALL LEAF NODE PAIRS IN THE GUIDE TREE MUST BE ASSOCIATED WITH AT LEAST 3 SEQUENCES!")
        print("\tTHE FOLLOWING LEAF NODE PAIRS DO NOT MEET THIS CRITERIA:\n")
        for pair in unsuitable_pairs:
            print(pair)
        print()
        comp_status = False

    else:
        print("\t[*] The Guide Tree, Imap, and Sequence Alignment are suitable for GDI calculations\n")
        comp_status = True

    return comp_status
