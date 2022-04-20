## DEPENDENCDIES
# HELPER FUNCTION DEPENDENCIES
from helper_functions import pretty
from helper_functions import pretty_Table
from helper_functions import read_MasterControl

# CHECKING FUNCTION DEPENDENCIES
from checking_functions import assert_Imap_Seq_compat
from checking_functions import check_Newick
from checking_functions import check_Numeric
from checking_functions import check_Tauprior
from checking_functions import check_Thetaprior
from checking_functions import check_GDI_params
from checking_functions import check_Finetune
from checking_functions import check_Threads
from checking_functions import check_BPP_mode
from checking_functions import check_Imap_Seq_compat
from checking_functions import check_Imap_filetype
from checking_functions import check_MSA_filetype
from checking_functions import check_Master_control_filetype
from checking_functions import check_Outfilename
from checking_functions import check_Print
from checking_functions import check_SandT_popsizes
from checking_functions import check_Imap_Tree_compat
from checking_functions import check_ValueIsFrom
from checking_functions import check_BPP_ctl_filetype
from checking_functions import check_nloci_MSA_compat
from checking_functions import check_Threads_nloci_compat
from checking_functions import check_SandT_Imap_MSA_compat

## DATA DEPENDENCIES
from data_dicts import master_Control_feedback
from data_dicts import BPP_Control_feedback
from data_dicts import col_print

## TYPE HINTS
from custom_types import BPP_control_dict, BPP_control_file, BPP_mode

## MASTER CONTROL CHECKING FUNCTION

# check the user supplied master control file for any misspecifications
def check_Master_Control(
        input_control_file: BPP_control_file
                        ):

    # begin by checking if the file can be loaded at all
    print(f"\n{col_print.BLUE}INITAL CHECK OF MASTER CONTROL FILE (MCF) PARAMETERS{col_print.RESETC}\n")
    
    # check that master control file is present and only contains calls to known parameters
    check_Master_control_filetype(input_control_file)
    

    ## SHALLOW CHECKING OF PARAMETERS
    # continue by checking all the parameters one-by-one
    param = read_MasterControl(input_control_file)
    par_check = {key:0 for key in param}

    # parameters for the pipeline
    par_check["file_align"]     = check_MSA_filetype(param["file_align"])
    par_check["file_imap"]      = check_Imap_filetype(param["file_imap"])
    par_check["tree_start"]     = check_Newick(param["tree_start"])
    par_check["tree_HM"]        = check_Newick(param["tree_HM"])
    par_check["ctl_file_phylo"] = check_BPP_ctl_filetype(param["ctl_file_phylo"])
    par_check["ctl_file_delim"] = check_BPP_ctl_filetype(param["ctl_file_delim"])
    par_check["ctl_file_HM"]    = check_BPP_ctl_filetype(param["ctl_file_HM"])

    # parameters for the merge decisions
    par_check["mode"]           = check_ValueIsFrom(param["mode"], ["merge", "split"])
    par_check["GDI_thresh"]     = check_Numeric(param["GDI_thresh"], "f", min = 0, max = 1)
    par_check["generations"]    = check_Numeric(param["generations"], "i", min = 1, max = 1000000)
    par_check["mutationrate"]   = check_Numeric(param["mutationrate"], "f", min = 0, max = 1)
    par_check["HM_decision"]    = check_GDI_params(param["HM_decision"], par_check)
    
    # parameters passed to BPP instances
    par_check["seed"]           = check_Numeric(param["seed"], "i", min = -10000)
    par_check['thetaprior']     = check_Thetaprior(param['thetaprior'])
    par_check['tauprior']       = check_Tauprior(param['tauprior'])
    par_check['finetune']       = check_Finetune(param['finetune'])
    par_check["sampfreq"]       = check_Numeric(param["sampfreq"], "i", min = 0, max = 100)
    par_check["nsample"]        = check_Numeric(param["nsample"], "i", min = 1000)
    par_check["burnin"]         = check_Numeric(param["burnin"], "i", min = 200)
    par_check['threads']        = check_Threads(param['threads'])
    

    ## PRINT RESULTS
    par_names = [*param]
    values = list(param.values())
    feedback = [master_Control_feedback[key][par_check[key]] for key in par_names]
    table = [par_names, values, feedback]
    colnames = ["PARAMETER", "VALUE", "FEEDBACK"]
    pretty_Table(table, colnames)


    ## MAKE FINAL DECISION TO PROCEED OR NOT
    # if no immediate errors are found, return which parameters are provided and which are empty
    print("SUMMARIZED RESULTS:")
    error_n = sum(i < 0 for i in list(par_check.values()))
    if error_n == 0:
        print(f"\n\t[*] No errors found during inital check of: '{input_control_file}'")
    
    # if erroneous parameters are found, halt execution immediately
    elif error_n > 0:
        print(f"\n{error_n} [X] ERROR(S) FOUND IN: '{input_control_file}'. PLEASE READ THE FEEDBACK, AND CONSULT THE MANUAL!")
        print("\t\t\n-- EXITING PROGRAM --")
        exit()


# check if the user supplied parameters to various BPP modes are correctly specified
def check_BPP_cfile (
        param:              BPP_control_dict, 
        sourcedict:         dict, 
        BPP_mode:           BPP_mode
                    ) ->    bool:

    mode_desc = {"A01":"starting phylogeny inference", "A11": "starting delimitation inference", "A00":"hierarchical method"}
    print(f"{col_print.BLUE}\nCHECKING USER SUPPLIED BPP {BPP_mode} PARAMETERS USED DURING {mode_desc[BPP_mode].upper()}{col_print.RESETC}\n")

    # prepare empty checking dict, the popsizes row is ignored because that is only an internal row of the pipeline
    par_check = {key:0 for key in param if key != "popsizes"} ## ADD REMOVAL OF UNSUPPLIED PARAMETERS

    ## INITITAL MISSPECIFICATION CHECKING
    # check file related parameters
    par_check['seqfile']        = check_MSA_filetype(param["seqfile"])
    par_check['Imapfile']       = check_Imap_filetype(param["Imapfile"])
    par_check['outfile']        = check_Outfilename(param["outfile"])
    par_check['mcmcfile']       = check_Outfilename(param["mcmcfile"])
    par_check["usedata"]        = check_ValueIsFrom(param["usedata"], ["1"])      
    
    # check if the BPP mode flags are correct for the intended case
    par_check["speciesdelimitation"], par_check["speciestree"] = check_BPP_mode(param["speciesdelimitation"], param["speciestree"], BPP_mode)
  
    # check model related parameters
    par_check["species&tree"]   = check_SandT_popsizes(param["species&tree"], param["popsizes"])
    par_check['newick']         = check_Newick(param['newick'])
    par_check['nloci']          = check_Numeric(param['nloci'], "i", min = 0)       
    par_check["cleandata"]      = check_ValueIsFrom(param["cleandata"], ["0","1"])       
    par_check['thetaprior']     = check_Thetaprior(param['thetaprior'])
    par_check['tauprior']       = check_Tauprior(param['tauprior'])
    par_check['finetune']       = check_Finetune(param['finetune'])
    
    # check MCMC related parameters
    par_check['seed']           = check_Numeric(param["seed"], "i", min = -10000)
    par_check['print']          = check_Print(param["print"])
    par_check["burnin"]         = check_Numeric(param["burnin"], "i", min = 200) 
    par_check["sampfreq"]       = check_Numeric(param["sampfreq"], "i", min = 0, max = 100)
    par_check["nsample"]        = check_Numeric(param["nsample"], "i", min = 1000) 
    par_check['threads']        = check_Threads(param['threads'])

    # check that the number of loci requested is not greater than the amount available in the MSA
    if par_check["nloci"] == 1 and par_check["seqfile"] == 1:
        par_check["nloci"] = check_nloci_MSA_compat(param["nloci"], param["alignmentfile"])
    # chcek that the number of threads requested does not exceed the number of loci in the MSA
    if par_check["nloci"] == 1 and par_check["threads"] == 1:
        par_check["threads"] = check_Threads_nloci_compat(param["threads"], param["nloci"])


    ## PRINT RESULTS OF INITIAL MISSPECIFICATION CHECKING
    par_names = [par_name for par_name in par_check if param[par_name] != "?"]
    source = [sourcedict[par_name] for par_name in par_names]
    values = [param[par_name] for par_name in par_names]
    feedback = [BPP_Control_feedback[key][par_check[key]] for key in par_names]
        # add a ghost row if the species&tree parameter was supplied. The ghost row displays the population size parameters
    if "species&tree" in par_names:
        insert_location = par_names.index("species&tree")+1
        par_names.insert(insert_location, "")
        source.insert(insert_location, "")
        values.insert(insert_location, param["popsizes"])
        feedback.insert(insert_location-1, "")
    
    table = [par_names, source, values, feedback]
    colnames = ["PARAMETER", "SOURCE", "VALUE", "FEEDBACK"]
    pretty_Table(table, colnames, width_limit = [2, 36])


    ## DEEP COMPATIBILITY CHECKING
        # only attempt deep compatibility checks if the required paramters are present
    if par_check["Imapfile"] and (par_check["newick"] or par_check["seqfile"]):
        print("\n\tDEEP COMPATIBILITY CHECKING\n")
        # check that the imap file and the MSA file refer to the same set of individuals
        par_check["imap_seq_compat"]         = check_Imap_Seq_compat(par_check, param["Imapfile"], param["seqfile"])
        # check that the imap file and the newick tree refer to the same set of populations
        par_check["imap_tree_compat"]        = check_Imap_Tree_compat(par_check, param["Imapfile"], param["newick"])
        # check that the population names and sizes in the species&tree row match the MSA and the Imap
        par_check["s_and_t_imap_msa_compat"] = check_SandT_Imap_MSA_compat(par_check, param["species&tree"], param["popsizes"], param["Imapfile"], param["seqfile"])


    ## MAKE FINAL DECISION TO PROCEED OR NOT
    print("\tSUMMARIZED RESULTS:")
    error_n = sum(i < 0 for i in list(par_check.values()))
    if error_n == 0:
        correctly_specified = True
        print(f"\n\t[*] No errors found in the user supplied parameters of BPP {BPP_mode}")
        
    
    # if erroneous parameters are found, halt execution immediately
    elif error_n > 0:
        correctly_specified = False
        print(f"\n{error_n} [X] ERROR(S) FOUND IN USER SUPPLIED PARAMETERS! PLEASE READ THE FEEDBACK, AND CONSULT THE MANUAL!")
    
    return correctly_specified

# check that the specified parameters for A01 and A11 are not in conflict
def check_A01_to_A11_compatibility  (
        A01_parameters:                     BPP_control_dict, 
        A11_parameters:                     BPP_control_dict,
                                    ) ->    bool:

    print(f"\n{col_print.BLUE}CHECKING IF A01 PARAMETERS AND OUTPUT WILL BE COMPATIBLE WITH THE A11 STAGE{col_print.RESETC}\n")
    
    # this is only the case if the Imap file is identical
    if A01_parameters["Imapfile"] == A11_parameters["Imapfile"]:
        compatible = True
        print("\t[*] The Imap for the A01 and A11 stages is identical, compatibility is guaranteed!")

    else:
        compatible = False
        print("\t[X] ERROR: IMAP FILE FOR A01 AND A11 STAGES IS NOT IDENTICAL!\n")
        print(f"\tA01 Imap: {A01_parameters['Imapfile']}")
        print(f"\tA11 Imap: {A11_parameters['Imapfile']}")
        print("\n\tPlease use the same Imap file for the two stages!\n")

    return compatible

# check that the specified parameters for A11 and A00 are not in conflict
def check_A11_to_A00_compatibility  (
        A11_parameters:                     BPP_control_dict, 
        A00_parameters:                     BPP_control_dict
                                    ) ->    bool:

    print(f"\n{col_print.BLUE}CHECKING IF A11 PARAMETERS AND OUTPUT WILL BE COMPATIBLE WITH THE A00 STAGE{col_print.RESETC}\n")

    A11_msa_checked = check_MSA_filetype(A11_parameters["seqfile"])
    if A11_msa_checked < 1:
        compatible = False
        if   A11_msa_checked == 0:
            print(f"\t[X] ERROR: NO ALIGNMENT FILE PROVIDED FOR THE A11 STAGE")
        elif A11_msa_checked == -1:
            print(f"\t[X] ERROR: NO FILE EXISTS AT REQUESTED LOCATION: {A11_parameters['seqfile']}")
        elif A11_msa_checked == -2:
            print(f"\t[X] ERROR: THE FILE {A11_parameters['seqfile']} IS NOT A VALID PHYLIP MSA")
    else:
        # the two stages can be compatible if the MSA files are identical
        if A11_parameters["seqfile"] == A00_parameters["seqfile"]:
            compatible = True
            print("\t[*] The MSA for the A11 and A00 stages is identical, compatibility is guaranteed!")
        elif assert_Imap_Seq_compat(A11_parameters["Imapfile"], A00_parameters["alignmentfile"]):
            compatible = True
            print("\t[*] The MSA for the A11 and A00 stages is different, but they reference the same individuals, so comaptibility is guaranteed!")
        else:
            compatible = False

    return compatible