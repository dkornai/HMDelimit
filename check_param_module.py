## DEPENDENCDIES

# HELPER FUNCTION DEPENDENCIES
from helper_functions import pretty, read_MasterControl
from helper_functions import prettyTable

# CHECKING FUNCTION DEPENDENCIES
from checking_functions import check_BPP_ctl_filetype, check_File, check_Imap_Seq_compat, check_Imap_filetype, check_MSA_filetype, check_Outfilename, check_Print, check_SandT_popsizes, check_Imap_Tree_compat, check_ValueIsFrom, check_speciesdelimitation
from checking_functions import check_Newick
from checking_functions import check_Numeric
from checking_functions import check_Tauprior
from checking_functions import check_Thetaprior
from checking_functions import check_GDI_params
from checking_functions import check_Finetune
from checking_functions import check_Threads

# BPP CONTROL FILE READING FUNCTIONS
from bpp_cfile_module import get_known_BPP_param

# DATA DEPENDENCIES
from data_dicts import master_Control_feedback


## MASTER CONTROL CHECKING FUNCTION

# check the user supplied master control file for any immediately notiable error
def check_Master_Control(input_control_file):
    # begin by checking if the file can be loaded at all
    print("\nINITAL CHECK OF MASTER CONTROL FILE (MCF) PARAMETERS:\n")
    try:
        param = read_MasterControl(input_control_file)
        print(f"MCF succesfully read from: {input_control_file}\n")
    except:
        print(f"ERROR: NO MASTER CONTROL FILE FOUND AT: {input_control_file}")
        exit()

    # continue by checking all the parameters one-by-one
    param_checked = {key:0 for key in param}

    # parameters for the pipeline
    param_checked["file_align"]     = check_File(param["file_align"])
    if param_checked["file_align"] == 1:
        param_checked["file_align"] = check_MSA_filetype(param["file_align"])
    
    param_checked["file_imap"]      = check_File(param["file_imap"])
    if param_checked["file_imap"]  == 1:
        param_checked["file_imap"]  = check_Imap_filetype(param["file_imap"])
    
    param_checked["tree_start"]     = check_Newick(param["tree_start"])
    param_checked["tree_HM"]        = check_Newick(param["tree_HM"])
    
    param_checked["ctl_file_phylo"] = check_File(param["ctl_file_phylo"])
    if param_checked["ctl_file_phylo"] == 1:
        param_checked["ctl_file_phylo"] = check_BPP_ctl_filetype(param["ctl_file_phylo"])
    
    param_checked["ctl_file_delim"] = check_File(param["ctl_file_delim"])
    if param_checked["ctl_file_delim"] == 1:
        param_checked["ctl_file_delim"] = check_BPP_ctl_filetype(param["ctl_file_delim"])
    
    param_checked["ctl_file_HM"]    = check_File(param["ctl_file_HM"])
    if param_checked["ctl_file_HM"]    == 1:
        param_checked["ctl_file_HM"]    = check_BPP_ctl_filetype(param["ctl_file_HM"])

    # parameters for the merge decisions
    param_checked["mode"]           = check_ValueIsFrom(param["mode"], ["merge", "split"])
    param_checked["GDI_thresh"]     = check_Numeric(param["GDI_thresh"], "f", min = 0, max = 1)
    param_checked["generations"]    = check_Numeric(param["generations"], "i", min = 1, max = 1000000)
    param_checked["mutationrate"]   = check_Numeric(param["mutationrate"], "f", min = 0, max = 1)
    param_checked["HM_decision"]    = check_GDI_params(param["HM_decision"], param_checked)
    
    # parameters passed to BPP instances
    param_checked["seed"]           = check_Numeric(param["seed"], "i", min = -10000)
    param_checked['thetaprior']     = check_Thetaprior(param['thetaprior'])
    param_checked['tauprior']       = check_Tauprior(param['tauprior'])
    param_checked['finetune']       = check_Finetune(param['finetune'])
    param_checked["sampfreq"]       = check_Numeric(param["sampfreq"], "i", min = 1, max = 1000000)
    param_checked["nsample"]        = check_Numeric(param["nsample"], "i", min = 1000)
    param_checked["burnin"]         = check_Numeric(param["burnin"], "i", min = 250)
    param_checked['threads']        = check_Threads(param['threads'])
    
    # printout of results
    parameters = [*param]
    values = list(param.values())
    feedback = []
    for key in param:
        feedback.append(master_Control_feedback[key][param_checked[key]])
    table = [parameters, values, feedback]
    colnames = ["PARAMETER", "VALUE", "FEEDBACK"]
    prettyTable(table, colnames)

    # if no immediate errors are found, return which parameters are provided and which are empty
    error_n = sum(i < 0 for i in list(param_checked.values()))
    
    if error_n == 0:
        print(f"\nNo errors found during inital check of: '{input_control_file}'")
        print("\n\t\t-- PROCEEDING TO NEXT STEP --")

        return param_checked
    
    # if erroneous parameters are found, halt execution immediately
    elif error_n > 0:
        print(f"\n{error_n} ERROR(S) FOUND IN: '{input_control_file}'. PLEASE READ THE FEEDBACK, AND CONSULT THE MANUAL!")
        print("\t\t\n-- EXITING PROGRAM --")
        exit()

def check_BPP_cfile(input_control_file, BPP_mode):
    mc_dict = read_MasterControl(input_control_file)
    param = get_known_BPP_param(mc_dict, BPP_mode, user_only= True)

    # prepare empty checking dict
        # the popsizes row is ignored because that is only an internal row of the pipeline
        # in BPP, that row is just associated with "species&tree"
    param_checked = {key:0 for key in param if key != "popsizes"}

    ## MISSPECIFICATION CHECKING

    # check seed
    param_checked['seed'] = check_Numeric(param["seed"], "i", min = -10000)
    
    # check if seq file exists, and check if seq file is actually a seq file
    param_checked['seqfile'] = check_File(param["seqfile"])
    if param_checked["seqfile"] == 1:
        param_checked["seqfile"] = check_MSA_filetype(param["seqfile"])

    # check if imap file exists, and check if the file is actually an imap file
    param_checked['Imapfile'] = check_File(param["Imapfile"])
    if param_checked["Imapfile"] == 1:
        param_checked["Imapfile"] = check_Imap_filetype(param["Imapfile"])
                          

    param_checked['outfile']        = check_Outfilename(param["outfile"])
    param_checked['mcmcfile']       = check_Outfilename(param["mcmcfile"])
    
    # check if the BPP mode flags are correct for the intended case
    if BPP_mode == "A01":
        param_checked["speciesdelimitation"] = check_ValueIsFrom(param["speciesdelimitation"], ["0"])
        param_checked["speciestree"]         = check_ValueIsFrom(param["speciestree"], ["1"])
    elif BPP_mode == "A11":
        param_checked["speciesdelimitation"] = check_speciesdelimitation(param["speciesdelimitation"])
        param_checked["speciestree"]         = check_ValueIsFrom(param["speciestree"], ["1"])
    elif BPP_mode == "A00":
        param_checked["speciesdelimitation"] = check_ValueIsFrom(param["speciesdelimitation"], ["0"])
        param_checked["speciestree"]         = check_ValueIsFrom(param["speciestree"], ["0"])              
    
    # check if the s&t and popsizes rows are correctly formatted, and contain no internal contradictions
    param_checked["species&tree"]   = check_SandT_popsizes(param["species&tree"], param["popsizes"])
                    
    # check other paramaters
    param_checked['newick']         = check_Newick(param['newick'])
    param_checked['nloci']          = check_Numeric(param['nloci'], "i", min = 0)       
    param_checked["usedata"]        = check_ValueIsFrom(param["usedata"], ["1"])                    
    param_checked["cleandata"]      = check_ValueIsFrom(param["cleandata"], ["0","1"])       
    param_checked['thetaprior']     = check_Thetaprior(param['thetaprior'])
    param_checked['tauprior']       = check_Tauprior(param['tauprior'])
    param_checked['finetune']       = check_Finetune(param['finetune'])
    param_checked['print']          = check_Print(param["print"])
    param_checked["burnin"]         = check_Numeric(param["burnin"], "i", min = 250) 
    param_checked["sampfreq"]       = check_Numeric(param["sampfreq"], "i", min = 1, max = 1000000)
    param_checked["nsample"]        = check_Numeric(param["nsample"], "i", min = 1000) 
    param_checked['threads']        = check_Threads(param['threads'])

    pretty(param_checked)

    ## DEEP COMPATIBILITY CHECKING
    
    # check if the Imap and the alignment are compatible
        # only attempt checking if these files are valid in the first place
    if param_checked["Imapfile"] == 1 and param_checked["seqfile"] == 1:
        comp_Imap_seq    = check_Imap_Seq_compat(param["Imapfile"], param["seqfile"])
    
    # check if the Imap and the supplied Newick tree are compatible
    if param_checked["Imapfile"] == 1 and param_checked["newick"] == 1:
        comp_Imap_newick = check_Imap_Tree_compat(param["Imapfile"], param["newick"])    


# print("@@@@@@\n")
# check_BPP_cfile("micro_test.txt", "A01")
# check_BPP_cfile("micro_test.txt", "A11")
# check_BPP_cfile("micro_test.txt", "A00")

print("@@@@@@\n")
check_BPP_cfile("badMC.txt", "A01")
check_BPP_cfile("badMC.txt", "A11")
check_BPP_cfile("badMC.txt", "A00")

# print("@@@@@@\n")
# check_BPP_cfile("Pcontrol.txt", "A01")
# check_BPP_cfile("Pcontrol.txt", "A11")
# check_BPP_cfile("Pcontrol.txt", "A00")

# print("@@@@@@\n")
# check_BPP_cfile("Pcontrol2.txt", "A01")
# check_BPP_cfile("Pcontrol2.txt", "A11")
# check_BPP_cfile("Pcontrol2.txt", "A00")