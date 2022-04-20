## DEPENDENCIES 

# PYTHON STANDARD DEPENDENCIES
import subprocess
import re
import os
import copy
import io
import warnings

# PYTHON LIBRARY DEPENDENCIES
import pandas as pd

from Bio import AlignIO

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

## DATA DEPENDENCIES
from data_dicts import empty_HM_parameters
from data_dicts import default_HM_parameters

## CORE HELPER FUNCTIONS
# sets the working directory to the source directory of the master control file, and also cuts the path.
def set_wd(mastercontrol_filepath):
    working_dir = os.path.dirname(mastercontrol_filepath)
    os.chdir(working_dir)
    print(f"\nWorking directory changed to: {working_dir}\n")

    def path_leaf(path):
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)

    return path_leaf(mastercontrol_filepath)

# reads a text file into an array of rows
def readLines(file_name):
    fileObj = open(file_name, "r")
    lines = fileObj.read().splitlines() 
    fileObj.close()
    
    return lines 

# return a flattened list from a list of lists
def flatten(t):
    return [item for sublist in t for item in sublist]

# pad a string to the required length using spaces
def padString(string, length):
    while len(string) < length:
            string = string + " "
    
    return string

# strip more than 1 trailing or leading whitespace
def stripall(input_string):
    result = input_string
    try:
        while result[0] == " " or result[-1] == " ":
            result = result.strip()
    except:
        pass

    return result

# strip out any rows with only whitespace characters
def removeEmptyRows(input_rows):  
    output = []
    for row in input_rows:
        # regex searches for non-whitespace characters
        if re.search("\S+", row):
            output.append(row)

    return output

# merges two dicts by updating the values of the first dict with several possible behaviours
def overwrite_dict(dict_1, dict_2,                                       # the two dicts
                   overwrite_to_unknown = False, ignore_known = False,):   # optional arguments to deviate from default behaviour
    dict_out = copy.deepcopy(dict_1)

    # overwrite any valiues in the old dict where the new dict is not "?"
    if overwrite_to_unknown == False and ignore_known == False:
        for key in dict_out:
            if key in dict_2 and dict_2[key] != "?":
                dict_out[key] = dict_2[key]
    
    # modifies default by overwriting for "?" values in the second dict
    elif overwrite_to_unknown == True and ignore_known == False:
        for key in dict_out:
            if key in dict_2:
                dict_out[key] = dict_2[key]
    
    # modifies the default behaviour, and does not overwrite for any parameters that are already known
    elif ignore_known == True and overwrite_to_unknown == False:
        for key in dict_out:
            if key in dict_2 and dict_1[key] == "?":
                dict_out[key] = dict_2[key]
    
    return dict_out 

# track the changes that occur when two dicts are merged, and return if a given parameter is from dict1 or dict2
def dict_source(dict_1, dict_2, dict_2_name, prev_sourcedict = None):
    if prev_sourcedict == None:
        sourcedict = {}
    else:
        sourcedict = copy.deepcopy(prev_sourcedict)
    
    for key in dict_2:
        if dict_2[key] != dict_1[key]:
            sourcedict[key] = dict_2_name
        else:
            if prev_sourcedict == None:
                sourcedict[key] = ""
                
    return sourcedict

## GENERAL FILE AND FOLDER I-O

# create a directory if it does not exist, and halt execution of the directory exists
def create_TargetDir(target_dir_name):
    try:
        os.mkdir(target_dir_name)
        print(f"Directory '{target_dir_name}' created") 
    
    except FileExistsError:
        print(f"ERROR: Directory '{target_dir_name}' already exists.")
        exit()


## PRINT HELPER FUNCTIONS

# print large control dicts in an easy-to-read way
def pretty(dict):
    longest_key_length = max(map(len, dict))
    
    for key in dict:
        printkey = str(key)
        while len(printkey) < longest_key_length:
            printkey = printkey + " "
        print(printkey,"=",dict[key])
    
    print()

# print a list of lists and associated column names in a readable form
def prettyTable(input_table, input_colnames, column_break = "  "):
    def customlen(string):
        string = str(string)
        string = string.split("\n")[0]
        
        return len(string)

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
def read_MCF_param(input_text, target_param):
    matches = str([match for match in input_text if target_param in match])
    try:
        result = matches.split("=")[1][:-2]
        try:
            result = result.split("#")[0]
        except:
            pass
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
def read_MasterControl(input_control_file):
    mc_lines = readLines(input_control_file)
    control_file_params = {# parameters for the pipeline
                           "file_align":      read_MCF_param(mc_lines, "alignment file"),
                           "file_imap":       read_MCF_param(mc_lines, "Imap file"), 
                           "tree_start":      read_MCF_param(mc_lines, "starting tree"),
                           "tree_HM":         read_MCF_param(mc_lines, "HM guide tree"),  
                           "ctl_file_phylo":  read_MCF_param(mc_lines, "BPP A01 starting phylogeny inference"),
                           "ctl_file_delim":  read_MCF_param(mc_lines, "BPP A11 starting delimitation"),           
                           "ctl_file_HM":     read_MCF_param(mc_lines, "BPP A00 HM parameter inference"),  
                           # parameters for the merge decisions
                           "mode":            read_MCF_param(mc_lines, "HM mode"),
                           "GDI_thresh":      read_MCF_param(mc_lines, "GDI threshold"),
                           "generations":     read_MCF_param(mc_lines, "HM generation threshold"),
                           "mutationrate":    read_MCF_param(mc_lines, "HM mutation"),
                           "HM_decision":     read_MCF_param(mc_lines, "HM decision parameters"),
                           # parameters passed to BPP instances
                           "seed":            read_MCF_param(mc_lines, "seed"),
                           "thetaprior":      read_MCF_param(mc_lines, "thetaprior"), 
                           "tauprior":        read_MCF_param(mc_lines, "tauprior"), 
                           "finetune":        read_MCF_param(mc_lines, "finetune"),
                           "sampfreq":        read_MCF_param(mc_lines, "sampfreq"), 
                           "nsample":         read_MCF_param(mc_lines, "nsample"),                            
                           "burnin":          read_MCF_param(mc_lines, "burnin"),
                           "threads":         read_MCF_param(mc_lines, "threads"),
                           }
    
    return control_file_params 

# given the values in the MCF, and default values, collect all necessary decision thresholds used in the HM stage
def get_HM_parameters(input_mc_dict):
    
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
def alignfile_to_MSA(align_file):
    align_raw = readLines(align_file)
    
    # remove all empty rows from the alignment, as this confuses biopython
    align_noempty = removeEmptyRows(align_raw)

    # remove superflous whitespaces from the numerical parameter rows, as this also confuses biopython
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
def Imap_to_List(imap_file):
    # ingest the imap file
    imap_rows = readLines(imap_file)
    
    # remove all empty rows from the IMAP, as this confuses the program needlessly
    imap_lines = removeEmptyRows(imap_rows)
    
    # create empty array to store results
    imap_indiv = ([i.split(None, 2)[0] for i in imap_lines])
    imap_pop = ([i.split(None, 2)[1] for i in imap_lines])
    
    return [imap_indiv, imap_pop]

# read the Imap text file or imap list to return a dict with key = individual ids, value = population assignments
def Imap_to_IndPop_Dict(imap, imap_is_list = False):
    # default behaviour, expects Imap files
    if imap_is_list == False:
        # ingest the imap file
        imap_rows = readLines(imap)
        # remove all empty rows from the IMAP, as this confuses the program needlessly
        imap_lines = removeEmptyRows(imap_rows)
         # create empty array to store results
        imap_indiv = [i.split(None, 2)[0] for i in imap_lines]
        imap_pop = [i.split(None, 2)[1] for i in imap_lines]
    
    # modified behaviour, handles Imap lists
    elif imap_is_list == True:
        imap_indiv = imap[0]
        imap_pop = imap[1]
    
    return dict(zip(imap_indiv, imap_pop))

# read the Imap text file, or an Imap list, to return a dict with key = pop names, value = IDs in that pop
def Imap_to_PopInd_Dict(imap, imap_is_list = False):
    # default behaviour, expects Imap files
    if imap_is_list == False:
        IDs = Imap_to_List(imap)[0]
        Population = Imap_to_List(imap)[1]
    # modified behaviour, handles Imap lists
    elif imap_is_list == True:
        IDs = imap[0]
        Population = imap[1]

    # categorize all IDs into populations
    popind_dict = {}
    for index, label in enumerate(Population):
        if label in popind_dict:
            popind_dict[label].append(IDs[index])
        else:
            popind_dict[label] = [IDs[index]]
    
    return popind_dict

# writes a nested list containing the columns of the Imap to a file
def list_To_Imap(imap_list, filename):
    output_file = '\n'.join('\t'.join(map(str,row)) for row in zip(imap_list[0],imap_list[1]))
    
    with open(filename, 'w') as f:
        f.write(output_file)

## BPP CFILE I-O FUNCTIONS

# read the parameters of a BPP control file into a dict
def bppcfile_to_dict(input_cfile):
    # import the control file to a dataframe
    df = pd.read_csv(input_cfile, sep = "=", header=None)
    
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
    
    # if a given parameter is present, but does not have a non-null value, return "?" for that parameter
    for index, row in df.iterrows():
        try: 
            if len(row["value"]) == 0:
                row["value"] = "?"
        except:
            row["value"] = "?"

    # make dataframe a dict
    dict = df.set_index('par').T.to_dict("records")[0]

    return(dict)

# write dict representing the BPP control file to disk
def dict_to_bppcfile(input_dict, cfile_name):
    # convert to pandas dataframe
    df = pd.DataFrame(list(input_dict.items()))
    
    # write dataframe to disk
    df.to_csv(cfile_name, sep = "=", header = False, index = False)
    with open(cfile_name, 'r+') as f:
        text = f.read()
        text = re.sub('popsizes=', '               ', text)
        text = re.sub('newick=', '               ', text)
        f.seek(0)
        f.truncate()
        f.write(text)


## BPP EXECUTABLE I-O FUNCTIONS

# run BPP with a given control file
def BPP_run(control_file):
    try:
        print("\nSTARTING BPP...\n")
        subprocess.run(["bpp", "--cfile", control_file])
        print("\n\t\t -- BPP RUN SUCCESSFUL --")
    
    except:
        print("ERROR: THE FILES SUPPLIED TO BPP CAUSE A CRASH")
        exit()

# extract the most probable species tree from the "outfile" produced by BPP A01 or BPP A11
def extract_Speciestree(control_file):
    try:
        # find the name of the output file
        BPP_outfile_name = bppcfile_to_dict(control_file)['outfile']
        # read in the output file
        input_file = readLines(BPP_outfile_name)
        
        # find the index of the row where the most probable trees are outputted, and add one to move to the trees
        rowindex_tree = [i for i, s in enumerate(input_file) if '(A)' in s][0]+1
        # extract the tree in newick format
        tree = re.search("\(.+\);" , input_file[rowindex_tree].split("  ")[-1]).group()
    
        return tree

    except:
        print("ERROR: TREE COULD NOT BE EXTRACTED FROM BPP RESULTS")
        exit() 

# extract the best supported species list produced by BPP A11
def extract_Pops(control_file):
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
        print("ERROR: POPULATIONS COULD NOT BE EXTRACTED FROM BPP RESULTS")
        exit() 
