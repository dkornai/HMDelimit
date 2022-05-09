'''
THIS MODULE CONTAINS HELPER FUNCTIONS THAT ARE REUSED IN OTHER MODULES.
MANY OF THE FUNCTIONS RELATE TO I-O OPERATIONS.
'''
## DEPENDENCIES 

# STANDARD LIBRARY DEPENDENCIES
from cgitb import text
import subprocess
import re
import os
import copy
import io
import warnings
import ntpath
import psutil

def kill(proc_pid):
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()

# EXTERNAL LIBRARY DEPENDENCIES
import pandas as pd

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

## DATA DEPENDENCIES
from data_dicts import empty_HM_parameters
from data_dicts import default_HM_parameters
from data_dicts import MCF_param_dict
from data_dicts import clprnt

## TYPING HINTS
from custom_types import file_path
from custom_types import Text_rows_list
from custom_types import Imap_file
from custom_types import Imap_list
from custom_types import Tree_newick
from custom_types import BPP_control_dict
from custom_types import BPP_control_file
from custom_types import Master_control_dict
from custom_types import Master_control_file
from custom_types import Phylip_MSA_file
from custom_types import Population_list
from custom_types import HM_decision_parameters
from custom_types import BPP_out_file

## CORE HELPER FUNCTIONS

# return a flattened list from a list of lists
def flatten (
        t:          list[list]
            ) ->    list:

    return [item for sublist in t for item in sublist]

# pad a string to the required length using spaces
def padString   (
        string:         str, 
        length:         int
                ) ->    str:

    while len(string) < length:
            string = string + " "
    
    return string

# limit the maximum length of a string
def string_limit(
        string:         str,
        limit:          int,
                ) ->    str:

    if len(string) > limit:
        string = f"{string[:limit-3]}..."
    
    return string

# strip more than 1 trailing or leading whitespace
def stripall(
        input_string:   str
            ) ->        str:

    result = input_string
    try:
        while result[0] == " " or result[-1] == " ":
            result = result.strip()
    except:
        pass

    return result

# reads a text file into an array of rows
def readLines   (
        file_name:      file_path
                ) ->    Text_rows_list:

    fileObj = open(file_name, "r")
    lines = fileObj.read().splitlines() 
    fileObj.close()
    
    return lines 

# strip out any rows with only whitespace characters
def remove_empty_rows   (
        input_rows:             Text_rows_list
                        ) ->    Text_rows_list:  

    output = [row for row in input_rows if re.search("\S+", row)]

    return output

# read a text file, and return a filtered version with all text after "#" and "*" removed
def read_filter_comments(
        input_file:             file_path
                        ) ->    Text_rows_list:

    lines = readLines(input_file)
    lines = [line.split("#")[0] for line in lines]
    lines = [line.split("*")[0] for line in lines]
    lines = remove_empty_rows(lines)

    return lines

# merges two dicts by updating the values of the first dict with several possible behaviours
def overwrite_dict  (
        dict_1:             dict, 
        dict_2:             dict,
        ow_to_unknown:      bool = False, # optional arguments to deviate from default behaviour
        ignore_known:       bool = False,
                    ) ->    dict:   

    dict_out = copy.deepcopy(dict_1)

    # overwrite any valiues in the old dict where the new dict is not "?"
    if ow_to_unknown == False and ignore_known == False:
        for key in dict_out:
            if key in dict_2 and dict_2[key] != "?":
                dict_out[key] = dict_2[key]
    
    # modifies default by overwriting for "?" values in the second dict
    elif ow_to_unknown == True and ignore_known == False:
        for key in dict_out:
            if key in dict_2:
                dict_out[key] = dict_2[key]
    
    # modifies the default behaviour, and does not overwrite for any parameters that are already known
    elif ignore_known == True and ow_to_unknown == False:
        for key in dict_out:
            if key in dict_2 and dict_1[key] == "?":
                dict_out[key] = dict_2[key]
    
    return dict_out 


## GENERAL FILE AND FOLDER I-O

# get just the filename part of a full file path
def path_filename   (
        path:               file_path
                    ):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

# create a directory if it does not exist, and halt execution of the directory exists
def create_TargetDir(
        target_dir_name:    file_path,
        user_message:       str = None
                    ):
    
    try:
        os.mkdir(target_dir_name)
        if user_message == None:
            print(f"Directory '{target_dir_name}' created")
        else:
            print(user_message)
    
    except FileExistsError:
        print(f"ERROR: Directory '{target_dir_name}' already exists.")
        exit()

## PRINT HELPER FUNCTIONS

# print large dicts in an easy-to-read way
def pretty  (
        dict:   dict
            ):
    
    longest_key_length = max(map(len, dict))
    
    for key in dict:
        printkey = str(key)
        while len(printkey) < longest_key_length:
            printkey = printkey + " "
        print("  ", printkey,"=",string_limit(str(dict[key]), 84))
    
    print()

# print a list of lists and associated column names in a readable form
def pretty_Table(   
        input_table:    list, 
        input_colnames: list, 
        column_break:   str = "  ", 
        width_limit:    list[int] = None
                ):

    def customlen(string): # custom function to ignore the part of the row after "\n" when calculating lengths
        return len(str(string).split("\n")[0])

    # limit the length of phrases in a column to a maximum if required
    if width_limit != None:
        target_col = width_limit[0]
        target_len = width_limit[1]

        for i in range(len(input_table[0])):
            item = input_table[target_col][i]
            if customlen(item) > target_len:
                input_table[target_col][i] = f"{item[:target_len-3]}..."

    # merge the column titles in as the first column elements
    for i in range(len(input_table)):
        input_table[i].insert(0, input_colnames[i])

    # find the width of the column by gettin the maximum length of the column elements
    col_width = []
    for column in input_table:
        col_width.append(max(map(customlen, column)))

    # print row by row
    for i in range(len(input_table[0])):
        toprint = ""
        for j in range(len(input_table)):
            toprint += (padString(input_table[j][i], col_width[j]) + column_break)
        
        print(toprint)    


## MASTER CONTROL FILE I-O FUNCTIONS

# this function is used to extract a named paramter from the control file, or return "?" if not found
def read_MCF_param  (
        input_text:         Text_rows_list, 
        target_param:       str
                    ) ->    str:

    matches = str([match for match in input_text if target_param in match])
    try:
        result = matches.split("=")[1][:-2]
        try:
            while result[0] == " " or result[-1] == " ":
                result = result.strip()
        except:
            if len(result) == 0:
                result = "?"        
    except:
        result = "?"

    return result

# extract all supplied parameters from the master control file to a dict
def read_MasterControl  (
        input_mc_file:          Master_control_file
                        ) ->    Master_control_dict:

    mc_lines = read_filter_comments(input_mc_file)

    # try to find lines that match the keyphrases of the parameters
    control_file_params = {param:read_MCF_param(mc_lines, MCF_param_dict[param]) for param in MCF_param_dict}
    
    return control_file_params 

# given the values in the MCF, and default values, collect all necessary decision thresholds used in the HM stage
def get_HM_parameters   (
        input_mc_dict:          Master_control_dict
                        ) ->    HM_decision_parameters:
    
    # collect empty and default values
    hm_par = overwrite_dict(empty_HM_parameters, default_HM_parameters)

    # overwrite with any information provided in the control file
    hm_par = overwrite_dict(hm_par, input_mc_dict)
    
    # if the user has not provided thresholds, place in the default ones
    thresh_dict = {"GDI_thresh":  "?",
                   "generations": "?"}
    if hm_par["mode"] == "merge":
        thresh_dict["GDI_thresh"] = default_HM_parameters["GDI_thresh_merge"]
        thresh_dict["generations"] = default_HM_parameters["max_gen"]
    elif hm_par["mode"] == "split":
        thresh_dict["GDI_thresh"] = default_HM_parameters["GDI_thresh_split"]
        thresh_dict["generations"] = default_HM_parameters["min_gen"]
    hm_par = overwrite_dict(hm_par, thresh_dict, ignore_known = True)

    # convert parameters to numbers to make comparisons viable
    hm_par["GDI_thresh"] = float(hm_par["GDI_thresh"])
    hm_par["generations"] = int(hm_par["generations"])
    if hm_par["mutationrate"] != "?":
        hm_par["mutationrate"] = float(hm_par["mutationrate"])

    return hm_par


## SEQUENCE ALIGNMENT I-O FUNCTIONS

# return a properly filtered BioPython MSA object when pointed to a valid alignment file
def alignfile_to_MSA(
        align_file:         Phylip_MSA_file
                    ) ->    list[MultipleSeqAlignment]:

    align_raw = readLines(align_file)
    
    # remove all empty rows, and superflous whitespaces from the numerical parameter rows, as this confuses biopython
    align_noempty = remove_empty_rows(align_raw)
    align_fixednum = []
    for line in align_noempty:
        # only applies to numeric lines, ignores others
        if re.search("[0-9]+[\s]+[0-9]+", line):
            line_fix = re.sub(" +", " ", line)
            line_fix = line_fix.strip()
            align_fixednum.append(line_fix)
        else:
            align_fixednum.append(line)
    
    # translate list to stringIO object that can be parsed by biopython
    align_str = "\n".join(align_fixednum)
    buf = io.StringIO(align_str)

    # read in the data from the buffered filtered version
    alignment_list = list(AlignIO.parse(buf,"phylip-relaxed")) #wrapped in list as AlignIO.parse returns generator
    
    return alignment_list


## IMAP I-O FUNCTIONS

# read the Imap text file to return a list with the individual ids and population assignments
def Imap_to_List(
        imap_file:      Imap_file
                ) ->    Imap_list:

    # ingest the imap file and remove all empty rows
    imap_rows = readLines(imap_file)
    imap_lines = remove_empty_rows(imap_rows)
    
    # create empty array to store results
    imap_indiv = ([i.split(None, 2)[0] for i in imap_lines])
    imap_pop = ([i.split(None, 2)[1] for i in imap_lines])
    
    return [imap_indiv, imap_pop]

# read the Imap text file or imap list to return a dict with key = individual ids, value = population assignments
def Imap_to_IndPop_Dict (
        imap, 
                        ):

    # modified behaviour, handles Imap lists
    if type(imap) == list:
        imap_indiv = imap[0]
        imap_pop = imap[1]
    # default behaviour, expects Imap files
    else:
        # ingest the imap file
        imap_rows = readLines(imap)
        # remove all empty rows from the IMAP, as this confuses the program needlessly
        imap_lines = remove_empty_rows(imap_rows)
         # create empty array to store results
        imap_indiv = [i.split(None, 2)[0] for i in imap_lines]
        imap_pop = [i.split(None, 2)[1] for i in imap_lines]
    
    return dict(zip(imap_indiv, imap_pop))

# read the Imap text file, or an Imap list, to return a dict with key = pop names, value = IDs in that pop
def Imap_to_PopInd_Dict (
        imap, 
                        ):

    # modified behaviour, handles Imap lists
    if type(imap) == list:
        IDs = imap[0]
        Population = imap[1]
    # default behaviour, expects Imap files
    else:
        IDs = Imap_to_List(imap)[0]
        Population = Imap_to_List(imap)[1]

    # categorize all IDs into populations
    popind_dict = {}
    for index, label in enumerate(Population):
        if label in popind_dict:
            popind_dict[label].append(IDs[index])
        else:
            popind_dict[label] = [IDs[index]]
    
    return popind_dict

# writes a nested list containing the columns of the Imap to a file
def list_To_Imap(
        imap_list:      Imap_list, 
        new_file_name:  file_path
                ):

    output_file = '\n'.join('\t'.join(map(str,row)) for row in zip(imap_list[0],imap_list[1]))
    
    with open(new_file_name, 'w') as f:
        f.write(output_file)


## BPP CFILE I-O FUNCTIONS

# read the parameters of a BPP control file into a dict
def bppcfile_to_dict(
        input_cfile:        BPP_control_file
                    ) ->    BPP_control_dict:

    # strip the comments to avoid potential confusion of downstream modules
    lines = read_filter_comments(input_cfile)
    buf = io.StringIO("\n".join(lines))

    # import the control file to a dataframe
    df = pd.read_csv(buf, sep = "=", header=None)
    
    # rename columns, convert every cell to text, and remove trailing and leading whitespaces
    df.rename(columns={0: "par", 1: "value"}, inplace=True)
    df = df.applymap(lambda x: stripall(x))

    # if present, reconfigure the rows under "species&tree" to move the required data into the value column
    try:
        st_row = df.index[df["par"] == "species&tree"].tolist()[0]
        df.loc[st_row+1, "value"] = df.loc[st_row+1, "par"]
        df.loc[st_row+1, "par"] = "popsizes"
        df.loc[st_row+2, "value"] = df.loc[st_row+2, "par"]
        df.loc[st_row+2, "par"] = "newick"
    except:
        pass

    # make dataframe a dict, and transform all types to string
    bpp_cdict = df.set_index('par').T.to_dict("records")[0]
    bpp_cdict = {str(param):str(bpp_cdict[param]) for param in bpp_cdict}

    # if a given parameter is present, but but has an empty value, return "?" for that parameter
    for param in bpp_cdict:
        if len(bpp_cdict[param]) == 0:
            bpp_cdict[param] = "?"

    return bpp_cdict

# write dict representing the BPP control file to disk
def dict_to_bppcfile(
        input_dict:         BPP_control_dict, 
        new_file_name:      file_path
                    ) ->    BPP_control_file:

    # convert to pandas dataframe
    df = pd.DataFrame(list(input_dict.items()))
    
    # write dataframe to disk
    df.to_csv(new_file_name, sep = "=", header = False, index = False)
    with open(new_file_name, 'r+') as f:
        text = f.read()
        text = re.sub('popsizes=', '               ', text)
        text = re.sub('newick=', '               ', text)
        f.seek(0)
        f.truncate()
        f.write(text)


## BPP EXECUTABLE I-O FUNCTIONS

# run BPP with a given control file
def BPP_run (
        control_file:   BPP_control_file
            ):

    try:
        print(f"{clprnt.GREEN}\nSTARTING BPP...\n")
        subprocess.run(["bpp", "--cfile", control_file])
        print(f"{clprnt.end}")
    
    except:
        print(f"{clprnt.end}\n[X] ERROR: UNEXPECTED EXIT FROM BPP")
        exit()

# run BPP with a given control file, and capture the stdout results
def BPP_run_capture (
        control_file:   BPP_control_file,
        proc_id
            ):

    process = subprocess.Popen(f"bpp --cfile {control_file}", shell = True, bufsize = 1,
                           stdout=subprocess.PIPE, stderr = subprocess.STDOUT,encoding='utf-8', errors = 'replace' ) 

    extime = ""
    while True:
        realtime_output = process.stdout.readline()
        if realtime_output == '' and process.poll() is not None:
            break
        if realtime_output:
            text_out = f"{realtime_output}"
            if "%" in text_out:
                percent = text_out.split()[0]
                if ":" in text_out.split()[-1]:
                    extime = f"time {text_out.split()[-1]}        "
                print("avg progress", percent, extime, end = '\r')
            if "Writing checkpoint file out" in text_out:
                print(text_out)
                kill(process.pid)
                

# resume a BPP run from a checkpoint, and capture the stdout results
def BPP_resume_capture (
        chkpoint_file,
        proc_id
            ):

    process = subprocess.Popen(f"bpp --resume {chkpoint_file}", shell = True, bufsize = 1,
                           stdout=subprocess.PIPE, stderr = subprocess.STDOUT,encoding='utf-8', errors = 'replace' ) 

    extime = ""
    while True:
        realtime_output = process.stdout.readline()
        if realtime_output == '' and process.poll() is not None:
            break
        if realtime_output:
            text_out = f"{realtime_output}"
            if "%" in text_out:
                percent = text_out.split()[0]
                if ":" in text_out.split()[-1]:
                    extime = f"time {text_out.split()[-1]}        "
                print("avg progress", percent, extime, end = '\r')
            # kill when checkpoint file is written
            if "Writing checkpoint file out" in text_out:
                print(text_out)
                kill(process.pid)
                

# get the summary of a BPP run with a given control file
def BPP_summary (
        control_file:   BPP_control_file
            ):

    process = subprocess.run(["bpp", "--summary", control_file], stdout=subprocess.PIPE, encoding='utf-8')
    out_text = process.stdout

    return out_text


# extract the most probable species tree from the "outfile" produced by BPP A01 or BPP A11
def extract_Speciestree (
        control_file:           BPP_control_file
                        ) ->    Tree_newick:

    try:
        # find the name of the output file
        BPP_outfile_name = bppcfile_to_dict(control_file)['outfile']
        # read in the output file
        input_file = readLines(BPP_outfile_name)

        # find the index of the row where the most probable trees are outputted, and add one to move to the trees
        rowindex_tree = [i for i, s in enumerate(input_file) if '(A)' in s][0]+1
        
        # extract the tree in newick format if a normal tree output was produced
        try:
            tree = re.search("\(.+\);" , input_file[rowindex_tree].split("  ")[-1]).group()
            
            return tree
        
        # handle edge case where BPP lumps all species into a single species
        except:
            tree = re.search("\(.+\)" , input_file[rowindex_tree]).group()
            tree += ";"
            print('\n>> BPP A11 OUTPUT INDICATES THAT ALL SUPPLIED POPULATIONS ARE A SINGLE SPECIES!')
            
            return tree

    except:
        print("\n[X] ERROR: NO TREE COULD NOT BE EXTRACTED FROM BPP RESULTS")
        exit() 

# extract the best supported species list produced by BPP A11
def extract_Pops(
        control_file:   BPP_control_file
                ) ->    Population_list:

    try:
        # find the name of the output file
        BPP_outfile_name = bppcfile_to_dict(control_file)['outfile']
        # read in the output file
        input_file = readLines(BPP_outfile_name)

        # find the index of the row where the most probable trees are outputted, and add one to move to the trees
        rowindex_tree = [i for i, s in enumerate(input_file) if '(A)' in s][0]+1
        
        # extract the tree in newick format
        speciesnames = input_file[rowindex_tree].split("  ")[-2][1:-1].split(" ")
        
        return speciesnames 

    except:
        print("[X] ERROR: POPULATIONS COULD NOT BE EXTRACTED FROM BPP RESULTS")
        exit() 

# get a dict of node names:tau values and node names:theta values from the outfile
def extract_Name_TauTheta_dict  (
        BPP_outfile:                    BPP_out_file
                                ) ->    tuple[dict[str, float], dict[str, float]]:

    lines = readLines(BPP_outfile)
    relevant_index = lines.index("List of nodes, taus and thetas:") # find the line where the node labels are listed
    lines = lines[relevant_index+2:]
    tau_dict = {line.split()[3]:float(line.split()[1]) for line in lines if float(line.split()[1]) != 0} #0 values excluded because they only occur when tau is not estimated
    theta_dict = {line.split()[3]:float(line.split()[2]) for line in lines if float(line.split()[2]) != -1} #-1 values excluded because they only occur when theta is not estimated (one seq for the population)

    return tau_dict, theta_dict