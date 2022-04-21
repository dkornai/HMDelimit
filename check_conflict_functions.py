'''
THESE FUNCTIONS EXIST TO CHECK FOR CONFLICTS BETWEEN FILES THAT ARE OF THE CORRECT TYPE.
FOR EXAMPLE, AN IMAP FILE AND ALIGNMENT FILE MIGHT NOT BE COMPATIBLE, BECAUSE INDIVIDUAL IDS
EXIST IN THE ALIGNMENT THAT DO NOT EXIST IN THE IMAP. 
'''
## DEPENDENCDIES
# STANDAR LIBRARY DEPENDENCIES
import warnings

# EXTERNAL LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

# HELPER FUNCTION DEPENDENCIES
from helper_functions import alignfile_to_MSA
from helper_functions import Imap_to_List
from helper_functions import Imap_to_PopInd_Dict

# ALIGNMENT AND IMAP SPECIFIC DEPENDENCIES
from align_imap_module import autoPopParam
from align_imap_module import count_Seq_Per_Pop



# check if an Imap file and a sequence alignment are mutually compatible
'''
This fucntion checks that all individual IDs in the alignment are mapped 
to a population in the IMAP. If this is the case, BPP can proceed successfullu
'''
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
        if param_check["Imapfile"] == 1 and param_check["seqfile"] == 1:
            comp_status = assert_Imap_Seq_compat(imapfile, alignmentfile)

        else: # if the files exist but have errors, they cannot be checked for compatibility
            print("\t[X] ERROR: MSA-IMAP COMPATIBILTIY CANNOT BE ASSESED DUE TO UPSTREAM MISSPECIFICATIONS\n") 
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
        if param_check["newick"] == 1 and param_check["Imapfile"] == 1:
            comp_status = assert_Imap_Tree_compat(imapfile, tree)
        
        else: # if the files exist but have errors, they cannot be checked for compatibility
            print("\t[X] ERROR: TREE-IMAP COMPATIBILTIY CANNOT BE ASSESED DUE TO UPSTREAM MISSPECIFICATIONS\n") 
            comp_status = -1
    else: # if both or one of the files are missing, dont assess the validity
        comp_status = 0
    
    return comp_status


# check if the species&tree row implies the same populations structure that the Imap and the seq alignment do.
'''
This implies that the number of sequences that are associated with a population in the data, 
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
                    print(f"\t population: {pop} | # of sequences in seqfile: {auto_popparam_dict[pop]} | # suggested by user: {user_popparam_dict[pop]}")
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
        if param_check["species&tree"] == 1 and param_check["seqfile"] == 1 and param_check["Imapfile"] == 1 and param_check['imap_seq_compat'] == 1:
            comp_status = assert_SandT_Imap_MSA_compat(s_and_t, popsizes, imapfile, alignmentfile)

        else: # if the files exist but have errors, they cannot be checked for compatibility
            print("\t[X] ERROR: SPECIES&TREE-IMAP-MSA COMPATIBILTIY CANNOT BE ASSESED DUE TO UPSTREAM MISSPECIFICATIONS\n") 
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