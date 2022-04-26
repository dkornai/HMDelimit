'''
THESE FUNCTIONS ACT AS A SAFEGUARD TO THE REMAINDER OF THE PIPELINE.
THEY CHECK IF ANY OF THE PARAMETERS SUPPLIED BY THE USER WILL CAUSE
A DOWNSTREAM CRASH IN THE SYSTEM. THIS IS IMPORTANT, AS A CRASH MAY
OCCUR AFTER WEEKS OF RUNTIME, WHEN THE PIPELINE SWITCHES FROM ONE 
MODE OF BPP TO THE NEXT
'''
## DEPENDENCDIES
# HELPER FUNCTION DEPENDENCIES
from helper_functions import pretty_Table
from helper_functions import read_MasterControl

# CHECKING HELPER DEPENDENCIES
from check_helper_functions import check_Newick, check_folders_do_not_exist
from check_helper_functions import check_Numeric
from check_helper_functions import check_Tauprior
from check_helper_functions import check_Thetaprior
from check_helper_functions import check_GDI_params
from check_helper_functions import check_Finetune
from check_helper_functions import check_Threads
from check_helper_functions import check_BPP_mode
from check_helper_functions import check_Imap_filetype
from check_helper_functions import check_MSA_filetype
from check_helper_functions import check_Master_control_filetype
from check_helper_functions import check_Outfilename
from check_helper_functions import check_Print
from check_helper_functions import check_SandT_popsizes
from check_helper_functions import check_ValueIsFrom
from check_helper_functions import check_BPP_ctl_filetype
from check_helper_functions import check_nloci_MSA_compat
from check_helper_functions import check_Threads_MSA_compat
from check_helper_functions import check_Threads_nloci_compat
from check_helper_functions import check_locusrate

# CONFLICT CHECKING DEPENDENCIES
from check_conflict_functions import check_Imap_Tree_compat
from check_conflict_functions import check_Imap_Seq_compat
from check_conflict_functions import assert_Imap_Seq_compat
from check_conflict_functions import check_SandT_Imap_MSA_compat

## DATA DEPENDENCIES
from data_dicts import master_Control_feedback
from data_dicts import BPP_Control_feedback
from data_dicts import clprnt

## TYPE HINTS
from custom_types import BPP_control_dict
from custom_types import Master_control_file
from custom_types import BPP_mode

## MASTER CONTROL CHECKING FUNCTION

# check the user supplied master control file for any misspecifications
'''
This function aims to check if the master control file is suitable for guiding the execution of the pipeline. 
The function:
    
    1) verifies that the master control file exists, and only contains named parameters that are known to the pipeline
    2) verifies that the parameter values are formatted according to expectations

Execution of the pipeline is halted immediately if errors are detected at any stage.
'''
def check_Master_Control(
        input_control_file: Master_control_file
                        ):

    print(f"\n{clprnt.BLUE}INITAL CHECK OF MASTER CONTROL FILE (MCF) PARAMETERS{clprnt.end}\n")
    
    ## CHECK MASTER CONTROL FILE INTEGRITY
    check_Master_control_filetype(input_control_file)
    param = read_MasterControl(input_control_file)

    # check that the target folders that the output will be written to do not exist
    check_folders_do_not_exist(input_control_file)

    par_check = {key:0 for key in param}

    ## SHALLOW CHECKING OF MISSPECIFICATION
    # parameters for the pipeline
    par_check["seqfile"]        = check_MSA_filetype(param["seqfile"])
    par_check["Imapfile"]       = check_Imap_filetype(param["Imapfile"])
    par_check["tree_start"]     = check_Newick(param["tree_start"])
    par_check["tree_HM"]        = check_Newick(param["tree_HM"])
    par_check["ctl_file_phylo"] = check_BPP_ctl_filetype(param["ctl_file_phylo"])
    par_check["ctl_file_delim"] = check_BPP_ctl_filetype(param["ctl_file_delim"])
    par_check["ctl_file_HM"]    = check_BPP_ctl_filetype(param["ctl_file_HM"])

    # parameters for the merge decisions
    par_check["mode"]           = check_ValueIsFrom(param["mode"], ["merge", "split"])
    par_check["GDI_thresh"]     = check_Numeric(param["GDI_thresh"], "0<x<1")
    par_check["generations"]    = check_Numeric(param["generations"], "100<=x<=1000000","i")
    par_check["mutationrate"]   = check_Numeric(param["mutationrate"], "0<x<1")
    par_check["HM_decision"]    = check_GDI_params(param["HM_decision"], par_check)
    
    # parameters passed to BPP instances
    par_check["seed"]           = check_Numeric(param["seed"], "-1<=x","i")
    par_check['thetaprior']     = check_Thetaprior(param['thetaprior'])
    par_check['tauprior']       = check_Tauprior(param['tauprior'])
    par_check['finetune']       = check_Finetune(param['finetune'])
    par_check["sampfreq"]       = check_Numeric(param["sampfreq"], "0<x<100", "i")
    par_check["nsample"]        = check_Numeric(param["nsample"], "1000<=x", "i")
    par_check["burnin"]         = check_Numeric(param["burnin"], "200<=x","i")
    par_check['threads']        = check_Threads(param['threads'])
    

    ## PRINT RESULTS
    par_names = [*param]
    values = list(param.values())
    feedback = [master_Control_feedback[key][par_check[key]] for key in par_names]
    pretty_Table(input_table = [par_names, values, feedback], 
                 input_colnames = ["PARAMETER", "VALUE", "FEEDBACK"])


    ## MAKE FINAL DECISION TO PROCEED OR NOT
    error_n = sum(i < 0 for i in list(par_check.values()))
    
    if   error_n == 0:
        print(f"\n[*] No errors found during inital check of {input_control_file}")
    
    elif error_n > 0:
        print(f"\n[X] {error_n} ERROR(S) FOUND IN: '{input_control_file}'. PLEASE READ THE FEEDBACK, AND CONSULT THE MANUAL!")
        exit()

# check if the user supplied parameters to various BPP modes are correctly specified
'''
This function performs a very deep check of all user supplied BPP parameters. The parameters are
collected upstream of this function using "get_user_BPP_param", which collects from the master 
control file, and any stage specific custom control file that is relevant. These parameters are:
    
    1) first checked for simple misspecifications (eg parameter is provided in the wrong format)
    2) checked for small internal conflicts (eg more threads are requested than loci)
    3) checked for deep compatibility (eg Imap contains populations not found in the tree)

If the collected parameters pass these checks, the pipeline is guaranteed to execute without
crashes of BPP. 
'''
def check_BPP_param (
        param:              BPP_control_dict, 
        sourcedict:         dict, 
        BPP_mode:           BPP_mode,
        after_A11:          bool = False,   # when this parameter is on, the lack of seq and imap files is not interpreted as an error
                    ) ->    bool:

    mode_desc = {"A01":"starting phylogeny inference", "A11": "starting delimitation inference", "A00":"hierarchical method"}
    print(f"{clprnt.BLUE}\nCHECKING USER SUPPLIED BPP {BPP_mode} PARAMETERS USED DURING {mode_desc[BPP_mode].upper()}{clprnt.end}\n")

    # prepare empty checking dict, the popsizes row is ignored because that is only an internal row of the pipeline
    par_check = {key:0 for key in param if key != "popsizes"}

    ## INITITAL MISSPECIFICATION CHECKING
    # check file related parameters
    par_check['seqfile']        = check_MSA_filetype(param["seqfile"])
    par_check['Imapfile']       = check_Imap_filetype(param["Imapfile"])
        # if after A11 is not true, a lack of Imap and MSA files is a fatal error
    if after_A11 == False: 
        if par_check["seqfile"]   == 0:
            par_check["seqfile"]   = -5 
        if par_check["Imapfile"]  == 0:
            par_check["Imapfile"]  = -5

    par_check['outfile']        = check_Outfilename(param["outfile"])
    par_check['mcmcfile']       = check_Outfilename(param["mcmcfile"])
    par_check["usedata"]        = check_ValueIsFrom(param["usedata"], ["1"])      
    
    # check if the BPP mode flags are correct for the intended case
    par_check["speciesdelimitation"], par_check["speciestree"] = check_BPP_mode(param["speciesdelimitation"], param["speciestree"], BPP_mode)
  
    # check model related parameters
    par_check["species&tree"]   = check_SandT_popsizes(param["species&tree"], param["popsizes"])
    par_check['newick']         = check_Newick(param['newick'])
    par_check['nloci']          = check_Numeric(param['nloci'], "0<x","i")       
    par_check["cleandata"]      = check_ValueIsFrom(param["cleandata"], ["0","1"])       
    par_check['thetaprior']     = check_Thetaprior(param['thetaprior'])
    par_check['tauprior']       = check_Tauprior(param['tauprior'])
    par_check['finetune']       = check_Finetune(param['finetune'])
    
    # check MCMC related parameters
    par_check['seed']           = check_Numeric(param["seed"], "-1<=x","i")
    par_check['print']          = check_Print(param["print"])
    par_check["burnin"]         = check_Numeric(param["burnin"], "200<=x","i")
    par_check["sampfreq"]       = check_Numeric(param["sampfreq"], "0<x<100", "i")
    par_check["nsample"]        = check_Numeric(param["nsample"], "1000<=x", "i")
    par_check['threads']        = check_Threads(param['threads'])
    par_check['locusrate']      = check_locusrate(param['locusrate'])

    # check that the number of loci requested is not greater than the amount available in the MSA
    if par_check["nloci"] and par_check["seqfile"]:
        par_check["nloci"] = check_nloci_MSA_compat(param["nloci"], param["seqfile"])
    # check that the number of threads requested does not exceed the number of loci in the MSA
    if par_check["threads"] and par_check["seqfile"]:
        par_check["threads"] = check_Threads_MSA_compat(param["threads"], param["seqfile"])
    # check that the number of threads requested does not exceed the number of loci requested by the user
    if par_check["threads"] and par_check["nloci"]:
        par_check["threads"] = check_Threads_nloci_compat(param["threads"], param["nloci"])

    ## PRINT RESULTS OF INITIAL MISSPECIFICATION CHECKING
    par_names = [par_name for par_name in par_check if param[par_name] != "?" or par_check[par_name] < 0]
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
    
    pretty_Table(input_table = [par_names, source, values, feedback], 
                 input_colnames = ["PARAMETER", "SOURCE", "VALUE", "FEEDBACK"], 
                 width_limit = [2, 36])

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
    error_n = sum(i < 0 for i in list(par_check.values()))
    if error_n == 0:
        correctly_specified = True
        print(f"\n[*] No errors found in the user supplied parameters of BPP {BPP_mode}")
        
    
    # if erroneous parameters are found, halt execution immediately
    elif error_n > 0:
        correctly_specified = False
        print(f"\n{error_n} [X] ERROR(S) FOUND IN USER SUPPLIED BPP PARAMETERS! PLEASE READ THE FEEDBACK, AND CONSULT THE MANUAL!")
    
    return correctly_specified


# check that the specified parameters for A01 and A11 are not in conflict
'''
This function guarantees that the species tree output of the A01 stage will be compatible 
with the following A11 stage. This can only be the case, if the Imap file for the two stages is identical.
In this case, all of the species named in the tree will be present in the A11 imap as well, and 
the A11 stage can proceed accordingly.
'''
def check_A01_to_A11_compatibility  (
        A01_parameters:                     BPP_control_dict, 
        A11_parameters:                     BPP_control_dict,
                                    ) ->    bool:

    print(f"\n{clprnt.BLUE}CHECKING IF A01 PARAMETERS AND OUTPUT WILL BE COMPATIBLE WITH THE A11 STAGE{clprnt.end}\n")
    
    # this is only the case if the Imap file is identical
    if A01_parameters["Imapfile"] == A11_parameters["Imapfile"]:
        compatible = True
        print("[*] The Imap for the A01 and A11 stages is identical, compatibility is guaranteed!")

    else:
        compatible = False
        print("\t[X] ERROR: IMAP FILE FOR A01 AND A11 STAGES IS NOT IDENTICAL!\n")
        print(f"\tA01 Imap: {A01_parameters['Imapfile']}")
        print(f"\tA11 Imap: {A11_parameters['Imapfile']}")
        print("\n\tPlease use the same Imap file for the two stages!\n")

    return compatible


# check that the specified parameters for A11 and A00 are not in conflict
'''
This function guarantees that the output of the A11 stage will be able to get
passed on to the A00 stage without any conflicts occuring. The A11 stage outputs
two objects, an Imap file, and a guide tree (corresponding to the most probable delimitation).
As such, A11 stage output can be passed to teh A00 stage if the seqfile for
the two stages only references identical individuals. That way, the seqfile, Imap, and
tree for the A00 stage will also be compatible. 
'''
def check_A11_to_A00_compatibility  (
        A11_parameters:                     BPP_control_dict, 
        A00_parameters:                     BPP_control_dict
                                    ) ->    bool:

    print(f"\n{clprnt.BLUE}CHECKING IF A11 PARAMETERS AND OUTPUT WILL BE COMPATIBLE WITH THE A00 STAGE{clprnt.end}\n")

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
            print("[*] The seqfile for the A11 and A00 stages is identical, compatibility is guaranteed!")
        elif assert_Imap_Seq_compat(A11_parameters["Imapfile"], A00_parameters["alignmentfile"]):
            compatible = True
            print("[*] The seqfile for the A11 and A00 stages is different, but they reference the same individuals, so comaptibility is guaranteed!")
        else:
            compatible = False

    return compatible