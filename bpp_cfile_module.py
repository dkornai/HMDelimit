## DEPENDENCDIES

# PYTHON STANDARD DEPENDENCIES
import copy
import random

# HELPER FUNCTION DEPENDENCIES
from helper_functions import dict_source
from helper_functions import overwrite_dict
from helper_functions import bppcfile_to_dict

# IMAP, MSA AND TREE SPECIFIC DEPENDENCIES
from align_imap_module import autoPopParam
from align_imap_module import autoPrior
from align_imap_module import autoStartingTree

# DATA DEPENDENCIES
from data_dicts import empty_BPP_cfile_dict
from data_dicts import default_BPP_cfile_dict
from data_dicts import default_BPP_A00_cfile_dict
from data_dicts import default_BPP_A01_cfile_dict
from data_dicts import default_BPP_A11_cfile_dict

# TYPE HINTING DEPENDENCIES
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

Each of the 20 BPP parameters that the pipeline requires are collected from one of 
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

2) BPP mode agnostic  parameters included in the MCF
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
'''
def get_known_BPP_param (
        input_mc_dict:          Master_control_dict, 
        mode:                   BPP_mode
                        ) ->    BPP_control_dict: 
    
    # start by combining the empty and default generic BPP control dicts
    BPP_cdict = overwrite_dict(empty_BPP_cfile_dict, default_BPP_cfile_dict)

    # get defaults specific to the provided mode
    if mode   == "A00":
        BPP_cdict = overwrite_dict(BPP_cdict, default_BPP_A00_cfile_dict)
    elif mode == "A01":
        BPP_cdict = overwrite_dict(BPP_cdict, default_BPP_A01_cfile_dict)
    elif mode == "A11":
        BPP_cdict = overwrite_dict(BPP_cdict, default_BPP_A11_cfile_dict)

    # find any generic BPP parameters available in the master control dict
    BPP_cdict = overwrite_dict(BPP_cdict, input_mc_dict)
    file_locations = {"seqfile" :input_mc_dict["file_align"],
                      "Imapfile":input_mc_dict["file_imap"]}
    BPP_cdict = overwrite_dict(BPP_cdict, file_locations)

    # if the master control dict points to a BPP control file, extract the available parameters
    if mode == "A01":
        stage = 'ctl_file_phylo'
    elif mode == "A11":
        stage = 'ctl_file_delim'
    elif mode   == "A00":
        stage = 'ctl_file_HM'

    try: 
        bpp_dict  = bppcfile_to_dict(input_mc_dict[stage])
        BPP_cdict = overwrite_dict(BPP_cdict, bpp_dict)
    except:
        pass
    
    # final overwrite with MCF trees if available
    if mode == "A11":
        if input_mc_dict["tree_start"] != "?": # harvest starting tree if available
             BPP_cdict["newick"] = input_mc_dict["tree_start"]
    elif mode   == "A00":
        if input_mc_dict["tree_HM"] != "?":    # harvest HM guide tree if available
             BPP_cdict["newick"] = input_mc_dict["tree_HM"]

    return BPP_cdict


# get the BPP paramters that are set by the user in any form
'''
This function only collects BPP paramters from the MCF, or from the stage specific BPP control file of the MCF.
The function also returns a source dict which tells the origin of the paramter (MCF or BPP cfile).
This is useful when checking compatibilities and misspecifications in the BPP command file,
as the user can be pointed to the source of the incompatible or badly specified parameter
'''
### FIX SOURCE IS INCORRECT !!
def get_user_BPP_param  (
        input_mc_dict:          Master_control_dict, 
        mode:                   BPP_mode,
        mask_A11:               bool = False # optional ability to mask any parameters that would be inherited from the A11 stage  
                        ) ->    tuple[BPP_control_dict, dict]:    
    
    # find any BPP parameters available in the master control dict
    BPP_cdict = overwrite_dict(empty_BPP_cfile_dict, input_mc_dict)
    file_locations = {"seqfile" :input_mc_dict["file_align"],
                      "Imapfile":input_mc_dict["file_imap"]}
    BPP_cdict = overwrite_dict(BPP_cdict, file_locations)

    # keep track of the origin of parameters
    sourcedict = dict_source(empty_BPP_cfile_dict, BPP_cdict, "MCF")

    # if the master control dict points to a BPP control file, extract the available parameters
    if mode == "A01":
        stage = 'ctl_file_phylo'
    elif mode == "A11":
        stage = 'ctl_file_delim'
        if input_mc_dict["tree_start"] != "?": # harvest starting tree if available
             BPP_cdict["newick"] = input_mc_dict["tree_start"]
    elif mode   == "A00":
        stage = 'ctl_file_HM'
        if input_mc_dict["tree_HM"] != "?":    # harvest HM guide tree if available
             BPP_cdict["newick"] = input_mc_dict["tree_HM"]
    
    try: 
        bpp_dict  = bppcfile_to_dict(input_mc_dict[stage])
        BPP_cdict = overwrite_dict(BPP_cdict, bpp_dict)
            # keep track of the origin of parameters
        sourcedict = dict_source(empty_BPP_cfile_dict, BPP_cdict, "BPP Cfile", sourcedict)
    except:
        pass

    # final overwrite with MCF trees if available
    tree_dict = {}
    if mode == "A11":
        if input_mc_dict["tree_start"] != "?": # harvest starting tree if available
            tree_dict = {"newick":input_mc_dict["tree_start"]}
            sourcedict["newick"] = "MCF"
            BPP_cdict = overwrite_dict(BPP_cdict, tree_dict)
    elif mode   == "A00":
        if input_mc_dict["tree_HM"] != "?":    # harvest HM guide tree if available
            tree_dict = {"newick":input_mc_dict["tree_HM"]}
            sourcedict["newick"] = "MCF"
            BPP_cdict = overwrite_dict(BPP_cdict, tree_dict)

    # if the option to mask the parameters that A00 inherits after running A11 is on, make the following parameters unknown
        # When A00 is run after A11, These parameters depend on the results of A11, which are unpredictable in advance.  
    if mask_A11 == True:
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
   -seeds using a randomly generated value
   -population structure data by reading the alignment and the imap
   -priors by implementing the method of Prof. Bruce Rannala

By generating these values, this function enables People unfamiliar 
with BPP and/or Bayesian phylogenetics can use the pipeline.
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