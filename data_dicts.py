# custom class from colored printing
class clprnt:
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    end = '\033[0m'


## MASTER CONTROL FILE RELATED DATA

# the pairs of MCF parameter names:keyphrases that are searched in the control file to read them
MCF_param_dict =    {
# parameters for the pipeline
"seqfile"       :"seqfile",
"Imapfile"      :"Imapfile", 
"tree_start"    :"starting tree",
"tree_HM"       :"HM guide tree",  
"ctl_file_phylo":"BPP A01 starting phylogeny inference",
"ctl_file_delim":"BPP A11 starting delimitation",           
"ctl_file_HM"   :"BPP A00 HM parameter inference",  
# parameters for the merge decisions
"mode"          :"HM mode",
"GDI_thresh"    :"GDI threshold",
"generations"   :"HM generation threshold",
"mutationrate"  :"HM mutation rate",
"HM_decision"   :"HM decision criteria",
# parameters passed to BPP instances
"seed"          :"seed",
"thetaprior"    :"thetaprior", 
"tauprior"      :"tauprior", 
"finetune"      :"finetune",
"sampfreq"      :"sampfreq", 
"nsample"       :"nsample",                            
"burnin"        :"burnin",
"threads"       :"threads",
"nloci"         :"nloci",
"locusrate"     :"locusrate",
"cleandata"     :"cleandata",
                    }

# the list of master control file parameters that can be specified from the command line
command_line_params = list(MCF_param_dict)
command_line_params.append("working_dir")       # working dir is an extra cmdline specific parameter not notmally used
command_line_params.remove("tree_start")        # tree data cannot be provided due to the fact that bash hates "()"
command_line_params.remove("tree_HM")

# contains the written feedback for the checker of the master control file
master_Control_feedback =   {
# parameters for the pipeline
"seqfile":         {-2:"[X] ERROR: THE SPECIFIED SEQFILE IS NOT A VALID MSA\n\n\t Please give the name of a correctly phylip formatted MSA file.\n",
                    -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                    0 :" ~  alignment file not specified",
                    1 :"[*] alignment file found",
                    },
"Imapfile":        {-2:"[X] ERROR: THE SPECIFIED FILE IS NOT A VALID IMAP\n\n\t Please give the name of a correctly formatted imap file.\n",
                    -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                    0 :" ~  Imap file not specified",
                    1 :"[*] Imap file found",
                    },
"tree_start":      {-3:"[X] ERROR: NAMING CONFLICT BETWEEN LEAF NODES AND INTERNAL NODES\n\n\t Leaf node names overlap with the internal node names. e.g. '((A, B), AB)' also has an internal node also called 'AB'! \n",
                    -2:"[X] ERROR: STARTING TREE CONTAINS POLYTOMIES\n\n\t Tree is correctly newick formatted, but BPP will not accept any trees with polytomies! \n",
                    -1:"[X] ERROR: STARTING TREE NOT A VALID NEWICK TREE\n\n\t Please supply a correctly formatted newick tree, or leave empty\n",
                    0 :" ~  starting tree not specified",
                    1 :"[*] starting tree correctly specified",
                    },
"tree_HM":         {-3:"[X] ERROR: NAMING CONFLICT BETWEEN LEAF NODES AND INTERNAL NODES\n\n\t Leaf node names overlap with the internal node names. e.g. '((A, B), AB)' also has an internal node also called 'AB'! \n",
                    -2:"[X] ERROR: HM GUIDE TREE CONTAINS POLYTOMIES\n\n\t Tree is correctly newick formatted, but BPP will not accept any trees with polytomies! \n",
                    -1:"[X] ERROR: HM GUIDE TREE NOT A VALID NEWICK TREE\n\n\t Please supply a correctly formatted newick tree, or leave empty\n",
                    0 :" ~  HM guide tree not specified",
                    1 :"[*] HM guide tree correctly specified",
                    },
"ctl_file_phylo":  {-2:"[X] ERROR: THE FILE CAN NOT BE INTERPRETED AS A BPP CONTROL FILE\n\n\t Please consult the BPP manual for advice on BPP control files, or leave empty\n",
                    -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                    0 :" ~  BPP A01 Starting phylogeny inference control file not specified",
                    1 :"[*] BPP A01 Starting phylogeny inference control file successfully found",
                    },
"ctl_file_delim":  {-2:"[X] ERROR: THE FILE CAN NOT BE INTERPRETED AS A BPP CONTROL FILE\n\n\t Please consult the BPP manual for advice on BPP control files, or leave empty\n",
                    -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                    0 :" ~  BPP A01 Starting delimitation control file not specified",
                    1 :"[*] BPP A11 Starting delimitation control file successfully found",
                    },   
"ctl_file_HM":     {-2:"[X] ERROR: THE FILE CAN NOT BE INTERPRETED AS A BPP CONTROL FILE\n\n\t Please consult the BPP manual for advice on BPP control files, or leave empty\n",
                    -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                    0 :" ~  BPP A00 HM Parameter Inference control file not specified",
                    1 :"[*] BPP A00 HM Parameter Inference control file successfully found",
                    }, 
# parameters for the merge decisions
"mode":            {-1:"[X] ERROR: HM MODE INCORRECTLY SPECIFIED\n\n\t Please specify as 'merge' or 'split', or leave empty\n",
                    0 :" ~  HM mode not specified",
                    1 :"[*] HM mode correctly specified"
                    },
"GDI_thresh":      {-1:"[X] ERROR: GDI THRESHOLD INCORRECTLY SPECIFIED\n\n\t Please use a numeric value between 0 and 1, or leave empty\n",
                    0 :" ~  GDI threshold not specified",
                    1 :"[*] GDI threshold correctly specified",
                    },
"generations":    { -1:"[X] ERROR: SPECIATION THRESHOLD IN GENERATIONS INCORRECTLY SPECIFIED\n\n\t Please use an integer value between 100 and 1000000, or leave empty\n",
                    0 :" ~  Speciation threshold in generations not specified",
                    1 :"[*] Speciation threshold in generations correctly specified",
                    },
"mutationrate" :   {-1:"[X] ERROR: SUBSTITUTIONS / SITE / YEAR INCORRECTLY SPECIFIED\n\n\t Please use a small numeric value between 0 and 1, or leave empty\n",
                    0 :" ~  Substitutions / site / generation not specified",
                    1 :"[*] Substitutions / site / generation corectly specified",
                    },
"HM_decision":     {-2:"[X] ERROR: HM DECISION CRITERIA INCORRECTLY SPECIFIED\n\n\t Decision criteria include age, but node age cannot be calculated without a specified mutation rate\n",
                    -1:"[X] ERROR: HM DECISION CRITERIA INCORRECTLY SPECIFIED\n\n\t Please use one of the decision criteria listed in the manual, or leave parameter empty\n",
                    0 :" ~  HM decision criteria not specified",
                    1 :"[*] HM decision criteria correctly specified",
                    },
# parameters passed to BPP instances
"seed":            {-1:"[X] ERROR: BPP SEED INCORRECTLY SPECIFIED\n\n\t Please use -1 or a positive integer value, or leave empty\n",
                    0 :" ~  seed not specified",
                    1 :"[*] BPP seed correctly specified",
                    },
"thetaprior":      {-1:"[X] ERROR: THETA PRIOR INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                    0 :" ~  BPP theta prior not specified",
                    1 :"[*] BPP theta prior correctly specified",
                    },
"tauprior":        {-1:"[X] ERROR: TAU PRIOR INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                    0 :" ~  BPP tau prior not specified",
                    1 :"[*] BPP tau prior correctly specified",
                    },
"finetune":        {-1:"[X] ERROR: FINETUNE INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                    0 :" ~  BPP finetune not specified",
                    1 :"[*] BPP finetune correctly specified",
                    },
"sampfreq":        {-1:"[X] ERROR: MCMC SAMPLING FREQUENCY INCORRECTLY SPECIFIED\n\n\t Please use an 100 => integer >= 1, or leave empty\n",
                    0 :" ~  MCMC sampling frequency not specified",
                    1 :"[*] MCMC sampling frequency specified",
                    },
"nsample":         {-1:"[X] ERROR: MCMC SAMPLE NUMBER INCORRECTLY SPECIFIED\n\n\t Please use an integer value >= 1000, or leave empty\n",
                    0 :" ~  MCMC sample number not specified",
                    1 :"[*] MCMC sample number correctly specified",
                    },                            
"burnin":          {-1:"[X] ERROR: MCMC BURN-IN SAMPLES INCORRECTLY SPECIFIED\n\n\t Please use an integer value >= 200, or leave empty\n",
                    0 :" ~  MCMC burn-in iterations not specified",
                    1 :"[*] MCMC burn-in iterations correctly specified",
                    },
"threads":         {-2:"[X] ERROR: MORE THREADS REQUESTED THAN AVAILABLE ON THE COMPUTER\n\n\t Please consult the BPP manual for how to set threading, set to an integer value < the number of CPU cores available\n",
                    -1:"[X] ERROR: BPP THREADING INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                    0 :" ~  BPP threading not specified",
                    1 :"[*] BPP threading correctly specified",
                    },
"nloci":           {-1:"[X] ERROR: NLOCI NOT A POSITIVE INTEGER\n\n\t Please set to an integer > 1, or leave empty\n",
                    0 :" ~  nloci not specified",
                    1 :"[*] nloci correctly specified",
                    },
"locusrate":       {-1:"[X] ERROR: LOCUSRATE ERRONEOUSLY SPECIFIED\n\n\t Please check the BPP manual for advice on the parameter, or leave empty\n",
                    0 :" ~  locusrate not specified",
                    1 :"[*] locusrate correctly specified",
                    },
"cleandata":       {-1:"[X] ERROR: CLEANDATA MUST BE '0' or '1'\n\n\t Please set according to prefernces about ambiguous sites and gaps in the alignment, or leave empty\n",
                    0 :" ~  cleandata not specified",
                    1 :"[*] cleandata correctly specified",
                    },
                        }

## DATA USED IN THE HIERARCHICAL METHOD SECTION
# the empty HM decision parameter dict 
empty_HM_parameters   = {
"mode":            "?",
"GDI_thresh":      "?",
"generations":     "?",
"mutationrate":    "?",
"HM_decision":     "?",
                        }

# the default values for HM decision parameters
default_HM_parameters = {
"mode":            "merge",  # the program will proceed in merge mode if the user does not specify otherwise
"GDI_thresh_merge":"0.2",    # popoulations are considered a single species if their GDI values are below this value
"GDI_thresh_split":"0.7",    # popoulations are considered two species if their GDI values are above this value
"max_gen":         "1000",   # populatioons are considered a single species if they separated less than this many generations ago
"min_gen":         "10000",  # populations are considered two species if they separated more than this many generations ago
"mutationrate":    "?",      # default mutation rate is unkown, and if none is provided, the program will not compute the age of the split
"HM_decision":     "all",    # by default, all measured parameters should be withinn their respective thresholds to accept a proposal
                        }

# all the possible decision criteria during the HM process
HM_decision_criteria  = [
"one_gdi",       # only one gdi is enough to accept a proposal
"both_gdis",     # both gdis need to be within their threshold to accept a proposal
"age",           # only age is enough to accept a proposal                        
"one_gdi_&_age", # age and one gdi needs to be within the threshold to accept a proposal
"any",           # a minimum of one of the parameters within the threshold to accept a proposal
"any_two",       # two of the parameters must be within their threshold to accept a proposal
"all",           # all parameters must be within their threshold to accept a proposal
"none",          # all proposals are accepted regardless of results, this is useful for exploring the GDI values under any possible delimitation
                        ]

# these descptions of all the possible decision criteria during the HM process will be printed to the user
HM_decision_criteria_description  = {
"one_gdi"        :"if at least one GDI value is within the threshold.",      
"both_gdis"      :"if both GDIs values are within the threshold.", 
"age"            :"only if the split occured more than the required number generations ago.",                     
"one_gdi_&_age"  :"if at least one GDI value is within the threshold, and the split occured more than the required number of generations ago",
"any"            :"if at least one measured parameter is within its respective threshold.",         
"any_two"        :"if at least of two measured parameters are within their respective thresholds.",    
"all"            :"if all measured parameters are within their thresholds.",
"none"           :"regardless of parameter values and thresholds.",
                                    }


## BPP CONTROL FILE RELATED DATA
# contains the list of valid BPP parameters
valid_BPP_param_names = [
'arch', 
'burnin', 
'checkpoint', 
'cleandata',
'diploid', 
'finetune', 
'gammaprior', 
'heredity', 
'Imapfile', 
'locusrate', 
'mcmcfile', 
'nloci', 
'nsample', 
'outfile', 
'print', 
'sampfreq', 
'scaling', 
'seed', 
'seqfile', 
'sequenceerror', 
'speciesdelimitation', 
'speciesmodelprior', 
'speciestree', 
'species&tree', 
'tauprior', 
'thetaprior', 
'usedata',
'threads',
                        ]

# contains the list of parameters that need to be present in a BPP control file
empty_BPP_cfile_dict =  {
'seed':                 '?',
'seqfile':              '?', 
'Imapfile':             '?', 
'outfile':              '?', 
'mcmcfile':             '?',
'speciesdelimitation':  '?',
'speciestree' :         '?',
'species&tree':         '?', 
'popsizes':             '?', 
'newick':               '?', 
'usedata':              '?', 
'nloci':                '?',
'locusrate':            '?', 
'cleandata':            '?', 
'thetaprior':           '?', 
'tauprior':             '?', 
'finetune':             '?', 
'print':                '?', 
'burnin':               '?', 
'sampfreq':             '?', 
'nsample':              '?', 
'threads':              '?'
                        }

# contains default values shared by all modes of BPP
default_BPP_param =    {
'seed':                 '?',
'seqfile':              '?', 
'Imapfile':             '?', 
'outfile':              '?', 
'mcmcfile':             '?',
'speciesdelimitation':  '?',
'speciestree' :         '?',
'species&tree':         '?', 
'popsizes':             '?', 
'newick':               '?', 
'usedata':              '?', 
'nloci':                '?',
'locusrate':            '0',
'cleandata':            '0', 
'thetaprior':           '?', 
'tauprior':             '?', 
'finetune':             '1: .01 .0001 .005 .0005 .2 .01 .01 .01', 
'print':                '1 0 0 0', 
'burnin':               '2000', 
'sampfreq':             '1', 
'nsample':              '20000', 
'threads':              '1 ',
                        }

# contains default values for parameters of the BPP A01 control file
default_BPP_A01_param =    {
'outfile':              'startphylo_out.txt', 
'mcmcfile':             'startphylo_mcmc.txt',
'speciesdelimitation':  '0',
'speciestree' :         '1',
'usedata':              '1', 
                                }

# contains default values for parameters of the A11 BPP control file
default_BPP_A11_param =    {
'outfile':              'startdelim_out.txt', 
'mcmcfile':             'startdelim_mcmc.txt',
'speciesdelimitation':  '1 1 2 1',
'speciestree' :         '1',
'usedata':              '1', 
                                }

# contains default values for parameters of the A00 BPP control file
default_BPP_A00_param =    {
'outfile':              'HM_out.txt', 
'mcmcfile':             'HM_mcmc.txt',
'speciesdelimitation':  '0',
'speciestree' :         '0',
'usedata':              '1', 
                                }

# final wrapper dict that enables mode flag based access to defaults
default_BPP_modespecific_param =    {
"A01":  default_BPP_A01_param,
"A11":  default_BPP_A11_param,
"A00":  default_BPP_A00_param,
                                    }  


# contains the written feedback for the checker of the master control file
BPP_Control_feedback = {
"seed":                {-1:"[X] ERROR: BPP SEED INCORRECTLY SPECIFIED\n\n\t Please use -1 or a positive integer value, or leave empty\n",
                        0 :"_",
                        1 :"[*] BPP seed correctly specified",
                        },
"seqfile":             {-5:"[X] ERROR: BPP CANNOT RUN WITHOUT A SPECIFIED SEQFILE\n\n\t Please specify a seqfile in the Master control file or stage specific control file\n",
                        -2:"[X] ERROR: THE SPECIFIED FILE IS NOT A VALID MSA\n\n\t Please give the name of a correctly phylip formatted MSA file.\n",
                        -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                        0 :"_",
                        1 :"[*] MSA file found",
                        },
"Imapfile":            {-5:"[X] ERROR: BPP CANNOT RUN WITHOUT A SPECIFIED IMAP\n\n\t Please specify an Imap in the Master control file or stage specific control file\n",
                        -2:"[X] ERROR: THE SPECIFIED FILE IS NOT A VALID IMAP\n\n\t Please give the name of a correctly formatted imap file.\n",
                        -1:"[X] ERROR: NO FILE OF ANY TYPE AT REQUESTED LOCATION\n\n\t Please give the name of a valid file, or leave empty\n",
                        0 :"_",
                        1 :"[*] Imap file found",
                        },
"outfile":             {-1:"[X] ERROR: OUTFILE NAME INCORRECTLY SPECIFIED\n\n\t Please set this parameter to a file name ending in .txt, or leave empty\n",
                        0 :"_",
                        1 :"[*] outfile name correctly specified",
                        },
"mcmcfile":            {-1:"[X] ERROR: MCMCFILE NAME INCORRECTLY SPECIFIED\n\n\t Please set this parameter to a file name ending in .txt, or leave empty\n",
                        0 :"_",
                        1 :"[*] mcmcfile name correctly specified",
                        },
"speciesdelimitation": {-4:"[X] ERROR: 'SPECIESDELIMITAION' FLAG SHOULD START WITH '1' FOR A11\n\n\t Please set check the BPP Manual for how to enable A11 mode, or leave empty\n",
                        -3:"[X] ERROR: 'SPECIESDELIMITAION' FLAG SHOULD BE '0' FOR A01\n\n\t Please set value to '0' in the A01 control file, or leave empty\n",
                        -2:"[X] ERROR: 'SPECIESDELIMITAION' FLAG SHOULD BE '0' FOR A00\n\n\t Please set value to '0' in the A00 control file, or leave empty\n",
                        -1:"[X] ERROR: 'SPECIESDELIMITAION' INCORRECTLY FORMATTED\n\n\t Please consult the BPP Manual for how to format correctly, or leave empty\n",
                        0 :"_",
                        1 :"[*] speciesdelimitation flag correctly set",
                        },
"speciestree":         {-4:"[X] ERROR: 'SPECIESTREE' FLAG SHOULD START WITH '1' FOR A11\n\n\t Please set value to '1' in the A11 control file, or leave empty\n",
                        -3:"[X] ERROR: 'SPECIESTREE' FLAG SHOULD BE '1' FOR A01\n\n\t Please set value to '1' in the A01 control file, or leave empty\n",
                        -2:"[X] ERROR: 'SPECIESTREE' FLAG SHOULD BE '0' FOR A00\n\n\t Please set value to '0' in the A00 control file, or leave empty\n",
                        -1:"[X] ERROR: 'SPECIESTREE' INCORRECTLY FORMATTED\n\n\t Please use '0' or '1' according to the intended mode, or leave empty\n",
                        0 :"_",
                        1 :"[*] speciestree flag corrrectly set",
                        },
"species&tree":        {-1:"[X] ERROR: 'SPECIES&TREE' INCORRECTLY FORMATTED\n\n\t Please consult the BPP Manual for how to format correctly, or leave empty\n",
                        0 :"_",
                        1 :"[*] species&tree correctly specified",
                        },                                     
"newick":              {-3:"[X] ERROR: NAMING CONFLICT BETWEEN LEAF NODES AND INTERNAL NODES\n\n\t Leaf node names overlap with the internal node names. e.g. '((A, B), AB)' also has an internal node also called 'AB'! \n",
                        -2:"[X] ERROR: TREE CONTAINS POLYTOMIES\n\n\t Tree is correctly newick formatted, but BPP will not accept any trees with polytomies! \n",
                        -1:"[X] ERROR: NOT A VALID NEWICK TREE\n\n\t Please supply a correctly formatted newick tree, or leave empty\n",
                        0 :"_",
                        1 :"[*] tree correctly specified",
                        },
"usedata":             {-1:"[X] ERROR: USEDATA MUST BE SET TO '1'\n\n\t Please set to 1, or leave empty\n",
                        0 :"_",
                        1 :"[*] usedata correctly specified",
                        },
"nloci":               {-2:"[X] ERROR: NLOCI LARGER THAN THE NUMBER OF LOCI IN THE MSA\n\n\t Please set 'nloci' to be at most the number of loci in the MSA, or leave empty\n",
                        -1:"[X] ERROR: NLOCI NOT A POSITIVE INTEGER\n\n\t Please set to an integer > 1, or leave empty\n",
                        0 :"_",
                        1 :"[*] nloci correctly specified",
                        },
"locusrate":           {-1:"[X] ERROR: LOCUSRATE ERRONEOUSLY SPECIFIED\n\n\t Please check the BPP manual for advice on the parameter, or leave empty\n",
                        0 :" ~ locusrate not specified",
                        1 :"[*] locusrate correctly specified",
                        },
"cleandata":           {-1:"[X] ERROR: CLEANDATA MUST BE '0' or '1'\n\n\t Please set according to prefernces about ambiguous sites and gaps in the alignment, or leave empty\n",
                        0 :"_",
                        1 :"[*] cleandata correctly specified",
                        },
"thetaprior":          {-1:"[X] ERROR: BPP THETA PRIOR INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                        0 :"_",
                        1 :"[*] BPP theta prior correctly specified",
                        },
"tauprior":            {-1:"[X] ERROR: BPP TAU PRIOR INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                        0 :"_",
                        1 :"[*] BPP tau prior correctly specified",
                        },
"finetune":             {-1:"[X] ERROR: BPP FINETUNE INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                        0 :"_",
                        1 :"[*] BPP finetune correctly specified",
                        },
"print":               {-1:"[X] ERROR: PRINT INCORRECTLY FORMATTED\n\n\t Please consult the BPP manual, or leave empty\n",
                        0 :"_",
                        1 :"print correctly specified",
                        },
"sampfreq":            {-1:"[X] ERROR: MCMC SAMPLING FREQUENCY INCORRECTLY SPECIFIED\n\n\t Please use an 100 => integer >= 1, or leave empty\n",
                        0 :"_",
                        1 :"[*] MCMC sampling frequency specified",
                        },
"nsample":             {-1:"[X] ERROR: MCMC SAMPLE NUMBER INCORRECTLY SPECIFIED\n\n\t Please use an integer value >= 1000, or leave empty\n",
                        0 :"_",
                        1 :"[*] MCMC sample number correctly specified",
                        },                            
"burnin":              {-1:"[X] ERROR: MCMC BURN-IN SAMPLES INCORRECTLY SPECIFIED\n\n\t Please use an integer value >= 200, or leave empty\n",
                        0 :"_",
                        1 :"[*] MCMC burn-in iterations correctly specified",
                        },
"threads":             {-4:"[X] ERROR: MORE THREADS REQUESTED THAN 'NLOCI'\n\n\t Please lower the number of threads requested\n",
                        -3:"[X] ERROR: MORE THREADS REQUESTED THAN LOCI IN THE MSA\n\n\t Please lower the number of threads requested\n",
                        -2:"[X] ERROR: MORE THREADS REQUESTED THAN AVAILABLE ON THE COMPUTER\n\n\t Please consult the BPP manual for how to set threading, set to an integer value < the number of CPU cores available\n",
                        -1:"[X] ERROR: BPP THREADING INCORRECRLY SPECIFIED\n\n\t Please consult the BPP manual, or leave empty\n",
                        0 :"_",
                        1 :"[*] BPP threading correctly specified",
                        },
                    }

## DATA USED IN THE PAIRWISE DISTANCE CALCULATION
# available characters:
avail_chars = set(["T", "C", "G", "A", "R", "Y", "W", "S", "M", "K", "H", "B", "D", "V"])

# this is the distance between two nucleotides, copied directly from "https://github.com/brannala/bpps/blob/master/src/PriorFunc.js"
distance_dict = {   
'AA':0.0,
'AC':1.0,
'AG':1.0,
'AT':1.0,
'CC':0.0,
'CG':1.0,
'CT':1.0,
'GG':0.0,
'GT':1.0,
'TT':0.0,
'AR':0.500,
'AY':1.0,
'AS':1.0,
'AW':0.500,
'AK':1.0,
'AM':0.500,
'AB':1.0,
'AD':0.333,
'AH':0.333,
'AV':0.333,
'CR':1.0,
'CY':0.500,
'CS':0.500,
'CW':1.0,
'CK':1.0,
'CM':0.500,
'BC':0.333,
'CD':1.0,
'CH':0.333,
'CV':0.333,
'GR':0.500,
'GY':1.0,
'GS':0.500,
'GW':1.0,
'GK':0.500,
'GM':1.0,
'BG':0.333,
'DG':0.333,
'GH':1.0,
'GV':0.333,
'RT':1.0,
'TY':0.500,
'ST':1.0,
'TW':0.500,
'KT':0.500,
'MT':1.0,
'BT':0.333,
'DT':0.333,
'HT':0.333,
'TV':1.0,
'RR':0.500,
'RY':1.0,
'RS':0.250,
'RW':0.250,
'KR':0.250,
'MR':0.250,
'BR':0.167,
'DR':0.333,
'HR':0.250,
'RV':0.333,
'YY':0.500,
'SY':0.250,
'WY':0.250,
'KY':0.250,
'MY':0.250,
'BY':0.333,
'DY':0.167,
'HY':0.333,
'VY':0.167,
'SS':0.500,
'SW':1.0,
'KS':0.250,
'MS':0.250,
'BS':0.333,
'DS':0.167,
'HS':0.167,
'SV':0.333,
'WW':0.500,
'KW':0.250,
'MW':0.250,
'BW':0.167,
'DW':0.333,
'HW':0.333,
'VW':0.167,
'KK':0.500,
'KM':1.0,
'BK':0.333,
'DK':0.333,
'HK':0.167,
'KV':0.167,
'MM':0.500,
'BM':0.167,
'DM':0.167,
'HM':0.333,
'MV':0.333,
'BB':0.333,
'BD':0.222,
'BH':0.222,
'BV':0.222,
'DD':0.333,
'DH':0.222,
'DV':0.222,
'HH':0.333,
'HV':0.222,
'VV':0.333,
'CA':1.0,
'GA':1.0,
'TA':1.0,
'GC':1.0,
'TC':1.0,
'TG':1.0,
'RA':0.500,
'YA':1.0,
'SA':1.0,
'WA':0.500,
'KA':1.0,
'MA':0.500,
'BA':1.0,
'DA':0.333,
'HA':0.333,
'VA':0.333,
'RC':1.0,
'YC':0.500,
'SC':0.500,
'WC':1.0,
'KC':1.0,
'MC':0.500,
'CB':0.333,
'DC':1.0,
'HC':0.333,
'VC':0.333,
'RG':0.500,
'YG':1.0,
'SG':0.500,
'WG':1.0,
'KG':0.500,
'MG':1.0,
'GB':0.333,
'GD':0.333,
'HG':1.0,
'VG':0.333,
'TR':1.0,
'YT':0.500,
'TS':1.0,
'WT':0.500,
'TK':0.500,
'TM':1.0,
'TB':0.333,
'TD':0.333,
'TH':0.333,
'VT':1.0,
'YR':1.0,
'SR':0.250,
'WR':0.250,
'RK':0.250,
'RM':0.250,
'RB':0.167,
'RD':0.333,
'RH':0.250,
'VR':0.333,
'YS':0.250,
'YW':0.250,
'YK':0.250,
'YM':0.250,
'YB':0.333,
'YD':0.167,
'YH':0.333,
'YV':0.167,
'WS':1.0,
'SK':0.250,
'SM':0.250,
'SB':0.333,
'SD':0.167,
'SH':0.167,
'VS':0.333,
'WK':0.250,
'WM':0.250,
'WB':0.167,
'WD':0.333,
'WH':0.333,
'WV':0.167,
'MK':1.0,
'KB':0.333,
'KD':0.333,
'KH':0.167,
'VK':0.167,
'MB':0.167,
'MD':0.167,
'MH':0.333,
'VM':0.333,
'DB':0.222,
'HB':0.222,
'VB':0.222,
'HD':0.222,
'VD':0.222,
'VH':0.222,
                }