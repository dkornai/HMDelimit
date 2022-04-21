'''
THIS MODULE CONTAINS FUNCTIONS USED FOR COLLECTING THE PARAMETERS OF
BPP CONTROL FILES FROM THE MASTER CONTROL FILE, BPP CONTROL FILES,
DEFAULT VALUES, AND AUTOMATICALLY GENERATED VALUES
'''
## DEPENDENCDIES

# STANDARD LIBRARY DEPENDENCIES
import copy
import random

# HELPER FUNCTION DEPENDENCIES
from helper_functions import overwrite_dict
from helper_functions import bppcfile_to_dict

# IMAP, MSA AND TREE SPECIFIC DEPENDENCIES
from align_imap_module import autoPopParam
from align_imap_module import autoPrior
from align_imap_module import autoStartingTree

# DATA DEPENDENCIES
from data_dicts import empty_BPP_cfile_dict
from data_dicts import default_BPP_param
from data_dicts import default_BPP_modespecific_param

## TYPE HINTS
from custom_types import BPP_mode
from custom_types import BPP_control_dict
from custom_types import Imap_list
from custom_types import Tree_newick
from custom_types import Master_control_dict


## SPECIALIZED FUNCTIONS

# extract BPP control file parameters from available data
'''
This function is extremely important for the functionality of the pipeline.
This is because it is responsible for collecting all of the BPP control file 
parameters which are required to successfully start a BPP instance. 

Each of the ~20 BPP parameters that the pipeline requires are collected from one of 
the 4 sources below, which are listed in their order of importance, 
from least to most important. Higher importance sources will overwrite 
parameters collected from lower importance sources. 

0) BPP mode agnostic default BPP parameters found in 'data_dicts'
    These include parameters such as "threads" or "nsample", which are
    shared by modes of BPP

1) BPP mode specific default parameters found in 'data_dicts'
    These include the "speciestree" or "species&tree" parameters which dictate
    the specific mode (A00, A01, A11)

  Together, these default values ensure that the user needs to specify the minimal
  number of parameters in the MCF

2) BPP mode agnostic parameters included in the MCF
    These include "seed", "thetaprior", "tauprior" and also the location of imap
    and alignment files, and even starting trees.

  These values are shared between all BPP instances, and enable the user to
  take fine grained control over the relevant parameters of BPP, without having to
  create seperate and redundant control files.

3) BPP mode specific parameters found in BPP control files that are specified
at the "ctl_file_phylo","ctl_file_delim","ctl_file_HM" parameters of the MCF.

  This enables users to take complete control over each stage of the process,
  and present different parameters to each mode of BPP. For example, this could
  be useful if the user wants to use more samples in the A00 stage than during A01

4) The starting tree, or guide tree from the master control file, if they are present.
    When finding paramters for the A11 mode, a non "?" value for the "tree_start" 
    parameter of the master control file will overwrite any tree in the dedicated A11 control file.
    When finding paramters for the A00 mode, a non "?" value for the "tree_HM" 
    parameter of the master control file will overwrite any tree in the dedicated A00 control file.

  This final behaviour enables the user to supply phylogenetic information to the pipeline
  without having to make dedicated BPP control files. This makes the pipeline significantly
  more user friendly
'''
def get_known_BPP_param (
        input_mc_dict:          Master_control_dict, 
        BPP_mode:               BPP_mode
                        ) ->    BPP_control_dict: 

    # 0), 1) start by combining the default generic parameters with the default mode specific ones.
    BPP_cdict = overwrite_dict(default_BPP_param, default_BPP_modespecific_param[BPP_mode])

    # 2) find any BPP parameters available in the master control dict
    BPP_cdict = overwrite_dict(BPP_cdict, input_mc_dict)

    
    # 3) if the master control dict specifies a custom control file for the appropriate stage
    stage_code = {"A01":'ctl_file_phylo', "A11":'ctl_file_delim',"A00":'ctl_file_HM'}
    if input_mc_dict[stage_code[BPP_mode]] != "?":
       # extract the available parameters then overwrite with those parameters
        user_BPP_cfile = bppcfile_to_dict(input_mc_dict[stage_code[BPP_mode]])
        BPP_cdict = overwrite_dict(BPP_cdict, user_BPP_cfile)
    
    # 4) final overwrite with MCF trees if available
    if   BPP_mode == "A11":
        if input_mc_dict["tree_start"] != "?": # harvest starting tree if available
             BPP_cdict["newick"] = input_mc_dict["tree_start"]
    elif BPP_mode == "A00":
        if input_mc_dict["tree_HM"] != "?":    # harvest HM guide tree if available
             BPP_cdict["newick"] = input_mc_dict["tree_HM"]

    return BPP_cdict



# get the BPP paramters that are set by the user in any form
'''
This function works similarly to "get_known_BPP_param", but only collects user supplied 
values, and ignores default ones. The user can provide values through the master control file,
or stage specific custom BPP control files. The collection of user specified parameters only
is useful in the pre-pipeline checking stage. Here, the user specified inputs need to be checked
to ensure correct formatting, and compatibility between the stages of the pipeline. As default 
values, and automatically generated values to the BPP control file are guaranteed to be correct,
it is unnecessary to check any of them. 

The function also returns a source dict which tells the origin of the paramter (MCF or BPP cfile).
This is useful when checking compatibilities and misspecifications in the BPP command file,
as the user can be pointed to the source of the incompatible or badly specified parameter.

The "after_A11" option is used when checking the user supplied parameters of the A00 stage,
when an A11 stage is also used before. In this case, the output from the A11 will always be
formatted by the pipeline to be compatible with the next A00 stage, so the checking of
those parameters becomes unnecessary. However, when runnin only the A00 stage, it is
necessary to check them.
'''
### FIX SOURCE IS INCORRECT !!
def get_user_BPP_param  (
        input_mc_dict:          Master_control_dict, 
        BPP_mode:               BPP_mode,
        after_A11:              bool = False # optional ability to mask any parameters that would be inherited from the A11 stage  
                        ) ->    tuple[BPP_control_dict, dict]:    
    
    # 0A) find any BPP parameters available in the master control dict
    BPP_cdict = overwrite_dict(empty_BPP_cfile_dict, input_mc_dict)
    # 0B) keep track of the origin of newly found parameters
    sourcedict = {param:(" - " if (BPP_cdict[param] == "?") else "MCF") for param in BPP_cdict }  

    # 1A) if the master control dict points to a BPP control file...
    stage_code = {"A01":'ctl_file_phylo', "A11":'ctl_file_delim',"A00":'ctl_file_HM'}
    if input_mc_dict[stage_code[BPP_mode]] != "?":
        # extract the available parameters then overwrite with those parameters
        user_BPP_cfile = bppcfile_to_dict(input_mc_dict[stage_code[BPP_mode]])
        BPP_cdict = overwrite_dict(BPP_cdict, user_BPP_cfile)
        # 1B) update the source dict
        for param in BPP_cdict:
            if param in user_BPP_cfile and user_BPP_cfile[param] != "?":
                sourcedict[param] = f"BPP {BPP_mode} Cfile"
        
    # 2AB) final overwrite with MCF trees if available
    if   BPP_mode == "A11":
        if input_mc_dict["tree_start"] != "?": # harvest starting tree if available
            BPP_cdict["newick"] = input_mc_dict["tree_start"]
            sourcedict["newick"] = "MCF"
    elif BPP_mode == "A00":
        if input_mc_dict["tree_HM"] != "?":    # harvest HM guide tree if available
            BPP_cdict["newick"] = input_mc_dict["tree_HM"]
            sourcedict["newick"] = "MCF"

    # 3) if the option to mask the parameters that A00 inherits after running A11 is on, make the following parameters unknown
        # When A00 is run after A11, These parameters depend on the results of A11, which are unpredictable in advance.  
    if after_A11 == True:
        BPP_cdict["Imapfile"] = "?"
        BPP_cdict["species&tree"] = "?"
        BPP_cdict["newick"] = "?"
        BPP_cdict["popsizes"] = "?"

    return BPP_cdict, sourcedict



# generate parameters with no available defaults or user supplied values
'''
This function is core to the pipeline. It is able to generate missing parameter values
that are data specific, and thus impossible to provide default values for.

It can generate missing:
   -'seed' using a randomly generated value
   -'species&tree' and 'nloci' lines using the actual paramters observed in the data
   -'tauprior' and 'thetaprior' by implementing the method of Prof. Bruce Rannala

By generating these values, this function enables the master control file, or the 
specialized BPP control files to be the absolute minimum length
'''
def generate_unkown_BPP_param   (
        input_control_dict:             BPP_control_dict
                                ) ->    BPP_control_dict:
    
    BPP_cdict = copy.deepcopy(input_control_dict)
    
    # generate seed if missing
    if BPP_cdict['seed'] == '?':
         BPP_cdict['seed'] = random.randint(1, 100000)
        
    # fill in population identity, size and loci numbers 
    if any(parameter == '?' for parameter in [BPP_cdict['species&tree'], BPP_cdict['popsizes'], BPP_cdict['nloci']]):
        popparam = autoPopParam(alignmentfile = BPP_cdict['seqfile'], 
                                imap          = BPP_cdict['Imapfile'])
        
        # overwrite s&t and popsizes together, as these always need to match the assumptions of eachother
        if BPP_cdict['species&tree'] == "?" or BPP_cdict['popsizes'] == "?": 
            BPP_cdict['species&tree'] = popparam['species&tree']
            BPP_cdict['popsizes']     = popparam['popsizes']
        
        # only overwrite loci number if unknown, as lower than maximum loci numbers can be used to check errors
        if BPP_cdict['nloci'] == "?":
            BPP_cdict['nloci']    = popparam['nloci']
    
    # generate priors if any are missing
    if any(parameter == '?' for parameter in [BPP_cdict['thetaprior'], BPP_cdict['tauprior']]):
        
        prior = autoPrior(alignmentfile = BPP_cdict['seqfile'], 
                          imapfile      = BPP_cdict['Imapfile'])
        
        # only overwrite priors one-by-one
        if BPP_cdict['thetaprior'] == '?':
            BPP_cdict['thetaprior'] = prior['thetaprior']
        if BPP_cdict['tauprior'] == '?':
            BPP_cdict['tauprior'] = prior['tauprior']

    return BPP_cdict



# generate the starting tree in the BPP A01 control file, if it is missing
'''
This function is only used at the very begenning of a BPP A01 analysis if no tree is
provided by the user (tree can be in MCF, or in any of the BPP control files)
'''
def generate_unknown_BPP_tree   (
        input_control_dict:             BPP_control_dict
                                ) ->    BPP_control_dict:

    BPP_cdict = copy.deepcopy(input_control_dict)

    if BPP_cdict['newick'] == '?':
        BPP_cdict['newick'] = autoStartingTree(alignmentfile = BPP_cdict['seqfile'], 
                                               imapfile      = BPP_cdict['Imapfile'])

    return BPP_cdict



# generate the parameters of the BPP A00 control file that change according to a new proposal
'''
This function aims to update the BPP control file to be compliant with the proposed 
changes in population structure (which are reperesented by the new tree and the new imap).
The population parameters are calculated using "autoPopParam".
'''

def proposal_compliant_BPP_param(
        input_control_dict:             BPP_control_dict, 
        prop_imap:                      Imap_list, 
        prop_imap_name:                 str, 
        prop_tree:                      Tree_newick
                                ) ->    BPP_control_dict:

    BPP_cdict = copy.deepcopy(input_control_dict)
   
    # generate the s&t and popsizes parameters corresponding to the proposed imap
    prop_pop_param = autoPopParam(imap          = prop_imap,
                                  alignmentfile = BPP_cdict["seqfile"])
    
    # create a dict which specifies the keys that will be overwritten
    BPP_proposed_param = {'species&tree': prop_pop_param["species&tree"], # species&tree values corresponding to the proposal
                          'popsizes'    : prop_pop_param["popsizes"],     # the popsizes corresponding to the proposal
                          'newick'      : prop_tree,                      # use the proposed tree topology
                          'Imapfile'    : prop_imap_name}                 # the proposed imap corresponding to the tree topology
    
    # overwrite with parameters specific to the proposal
    BPP_cdict = overwrite_dict(BPP_cdict, BPP_proposed_param)

    return BPP_cdict