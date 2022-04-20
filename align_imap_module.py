'''
THIS MODULE CONTAINS SPECIALIZED FUNCTIONS FOR INTERPRETING, 
AND ANALYSING THE CONTENTS OF SEQUENCE ALIGNMENTS AND IMAP FILES.
'''

## DEPENDENCDIES

# PYTHON STANDARD DEPENDENCIES
import re
import copy
import random
import warnings
from io import StringIO
from collections import Counter
from itertools import combinations

# PYTHON LIBRARY DEPENDENCIES
import numpy as np

import Bio
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.Consensus import majority_consensus
from Bio import Phylo
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

# HELPER FUNCTION DEPENDENCIES
from helper_functions import Imap_to_PopInd_Dict
from helper_functions import Imap_to_IndPop_Dict
from helper_functions import alignfile_to_MSA

# DATA DEPENDENCIES
from data_dicts import distance_dict
from data_dicts import avail_chars

# TYPE HINTING DEPENDENCIES
from custom_types import Newick_tree
from custom_types import Imap_file

## SPECIALIZED HELPER FUNCTIONS

# calculate the pairwise distance between two sequences at shared known characters
def pairwise_dist(seq_1:str, seq_2:str) -> float:
    # collect sites where both sequences have comparable IUPAC codes
    seq_1_correct = set([i for i, x in enumerate(seq_1) if x in avail_chars])
    seq_2_correct = set([i for i, x in enumerate(seq_2) if x in avail_chars])
    overlap = list(seq_1_correct.intersection(seq_2_correct))
    alignlist = [f"{seq_1[i]}{seq_2[i]}" for i in overlap]
    
    # at each site, use the lookup table to measure the distance
    dist_persite = [distance_dict[pair] for pair in alignlist]
    
    # get the overall average
    avg_dist = (sum(dist_persite)/len(dist_persite))

    return np.round(avg_dist, decimals = 4)


# return a list of all paiwise distances in an alignment
def distanceList(input_MSA) -> list[float]:
    # isolate only the sequence strings
    seqlist = [str(sequence.seq) for sequence in input_MSA]
    
    # procude the list corresponding to which two lists will be compared in which order
        # the convoluted order is implemented to match up with the BioPython "DistanceMatrix"
    seq_com = sorted(list(combinations(list(range(len(seqlist))), 2)), key=lambda x: x[1])

    # go through the list of combinations, and measure the distance for each
    dist_list = [0]*len(seq_com)
    for i in range(len(seq_com)):
        seq1 = seqlist[seq_com[i][0]]
        seq2 = seqlist[seq_com[i][1]]
        dist_list[i] = pairwise_dist(seq1, seq2)

    return dist_list


# produce a BioPython compatible "DistanceMatrix" by wrapping the "distanceList" function

    # This distance matrix is used in "autoStartingTree" to generate trees using upgma

def getDistanceMatrix(input_MSA) -> DistanceMatrix:
    dist_list = distanceList(input_MSA)

    name_list = [str(seq.id) for seq in input_MSA]

    # format matrix to comply with BioPython by adding 0s, and getting values in the correct order
    matrix = []
    for x in range(len(input_MSA)):
        inrow = [0]
        for _ in range(x):
            inrow.insert(-1, dist_list.pop(0))
        matrix.append(inrow)
    
    return DistanceMatrix(names=name_list, matrix=matrix)


# return a dict containing the maximum number of sequences at any loci for each population

    # This function is used in "autoPopParam" to generate the descriptive lines for the BPP
    # control file. The function is also used in "check_GuideTree_Imap_compat" to ensure
    # that GDI values can be calculated. 

def count_Seq_Per_Pop(input_popind_dict, input_MSA_list):
    # create empty dict to hold results
    maxcounts = dict.fromkeys(input_popind_dict, 0)
    
    for curr_locus in input_MSA_list:
        # filter out a list of individual ids at a given locus
        curr_id_list = [seq.id for seq in curr_locus]
        curr_id_list = [re.sub('[a-zA-Z_0-9.]+\^', '', id) for id in curr_id_list]
        # replace individual ids with their population code
        curr_pop_list = [[k for k, v in input_popind_dict.items() if id in v][0] for id in curr_id_list]
        # count the number of sequences associated with each population, and keep track of the highest value
        counts = dict(Counter(curr_pop_list))
        for key in counts:
            if counts[key] > maxcounts[key]:
                maxcounts[key] = counts[key]
    
    return maxcounts


## SPECIALIZED FUNCTIONS

# generate the lines of the control file corresponding to population numbers and sizes

    # This function reads the alignment and the IMAP, and extracts the population labels, 
    # and the number of sequences associated with that population in the alignment.

def autoPopParam(imap, alignmentfile, imap_is_list = False):
    popind_dict = Imap_to_PopInd_Dict(imap, imap_is_list)   
    alignment = alignfile_to_MSA(alignmentfile)
    
    # row describing the number and name of populations
    n_pops = str(len(popind_dict))
    pop_names = str(popind_dict.keys())[11:-2]
    pop_names = pop_names.replace(",",""); pop_names = pop_names.replace("'","")
    
    # row describing the ,maximum number of of sequences/loci for each population 
    maxcounts = count_Seq_Per_Pop(popind_dict, alignment)
    n_in_pop = str([maxcounts[key] for key in maxcounts])[1:-1].replace(","," ")

    # row describing the number of loci
    nloci = str(len(alignment))

    # final output, formatted to comply with BPP control dict standards
    rows = {"species&tree": f"{n_pops} {pop_names}", 
            "popsizes"    : n_in_pop, 
            "nloci"       : nloci}

    return rows

# automatically generates priors (using the method of Prof. Bruce Rannala)

    # Theta is calculated to give a mean which is the average of within population average pairwise distances
    # tau is calculated to give a mean that is the largest pairwise distance seen in the alignment. 
    # The output of this function is used in BPP control files
 
def autoPrior(imapfile:Imap_file, alignmentfile):
    alignment = alignfile_to_MSA(alignmentfile)
    
    indpop_dict = Imap_to_IndPop_Dict(imapfile)
    populations = list(set(indpop_dict.values()))

    ## THETA CALCULATION
    # calculation of average-average pairwise distances
    dist_pop = []
        # iterate through all populations
    for population in populations:
        per_locus_dist = []
        per_locus_len = []

        for locus in alignment:
            # add the sequences belonging to the current population to a temp aligment
            temp_aligment = Bio.Align.MultipleSeqAlignment([])
            for sequence in locus:
                id = sequence.id
                id = id.split("^")[1]
                if indpop_dict[id] == population:
                    temp_aligment.append(sequence)
                # append the length of sequences in the temp alignment
            length = temp_aligment.get_alignment_length()
            
            # the distance can only be calculated for more than 2 sequences
            if len(temp_aligment) >= 2:
                
                # get avergage distance within the temp alignment
                dist_list = distanceList(temp_aligment)
                avg_dist = np.average(dist_list)

                # append to final list
                per_locus_len.append(length)
                per_locus_dist.append(avg_dist)
        
        # calculate the locus length weigthed within population average
        if len(per_locus_dist) > 0:
            pop_avg = np.average(per_locus_dist, weights = per_locus_len)
            dist_pop.append(pop_avg)
        
    # final theta calculation    
    D = np.average(dist_pop)
    theta_alpha = 3
    theta_beta = np.round(2*D, decimals = 4)
    
    ## TAU CALCULATION
    # calculation of the maximum pairwise distance at each locus
    max_dist = []
    for locus in alignment:
        dist_list = distanceList(locus)
        maxval = np.max(dist_list)
        max_dist.append(maxval)

    M = np.max(max_dist)
    tau_alpha = 3
    tau_beta = np.round(2*M, decimals = 4)

    # formatting of results to comply with BPP control file standards
    priors = {"thetaprior": f"{theta_alpha} {theta_beta} e",
              "tauprior"  : f"{tau_alpha} {tau_beta}",}

    return priors


# generate a starting tree using distance methods

    # This function is only called if the program needs to run BPP A01 to generate a starting tree.
    # The tree output is not very correct, but offers a better starting point that a random tree.
    # This way, less computational resources are wasted during A01.

def autoStartingTree(imapfile:Imap_file, alignmentfile) -> Newick_tree:
    # associate all individuals with populations
    indpop_dict = Imap_to_IndPop_Dict(imapfile)

    # associate all populations with individuals
    popind_dict = Imap_to_PopInd_Dict(imapfile)
    n_pop = len(popind_dict)
    
    # read in the MSA to a list containing the alignment at each locus
    alignment_list = alignfile_to_MSA(alignmentfile)

    ## GENERATE A LIST OF TREES FOR EACH LOCI BY RANDOMLY SAMPLING ONE SEQUENCE FROM EACH POPULATION
    tree_list = []
    # begin to iterate through the loci
    for locus in alignment_list:
        # get the ids of all the sequences at the locus
        seq_ids_at_locus = [str(record.id).split("^")[1] for record in locus]
        
        # check if all the population are present at the locus, and ignore the locus if not
        pops_at_locus = list(set([str(indpop_dict[indiv_id]) for indiv_id in seq_ids_at_locus]))
        if len(pops_at_locus) < n_pop:
            continue

        # recreate the two-way assocations specific to this locus
        popmap_at_locus = {pop:[seq_id for seq_id in popind_dict[pop] if seq_id in seq_ids_at_locus] for pop in pops_at_locus}
        pop_dict_at_locus = {seq_id: indpop_dict[seq_id] for seq_id in seq_ids_at_locus }

        # select one random individual from each population
        selected_ids = [random.choice(popmap_at_locus[key]) for key in popmap_at_locus]

        # create the alignment containing a sequence from the random individual
        temp_align = Bio.Align.MultipleSeqAlignment([])
        for seqObj in locus:
            indiv_id = str(seqObj.id).split("^")[1]
            
            if indiv_id in selected_ids: 
                selected_ids.remove(indiv_id)
                row_obj = copy.deepcopy(seqObj)
                row_obj.id = pop_dict_at_locus[indiv_id]
                temp_align.extend([row_obj])

        # infer tree using custom distance methods
        dm = getDistanceMatrix(temp_align)
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)
        tree_list.append(tree)

    ## USE A CONSESUS APPROACH WITH THE GENERATED TREES TO GENERATE OUTPUT
    majority_tree = majority_consensus(tree_list)
        # format to newick, and resolve any polytomies
    treeIO = StringIO()
    Phylo.write([majority_tree], treeIO, "newick")
    t = Tree(treeIO.getvalue())
    t.resolve_polytomy(recursive=True)
    t = str(t.write(format=9))
    
    return t
