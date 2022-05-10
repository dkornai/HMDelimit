'''
THIS MODULE CONTAINS THE FUNCTION TO CHECK IF THE 
PHYLOGENETIC SIGNAL FOR THE TREE TOPOLOGY IS SUFFICIENTLY STRONG.
'''


## DEPENDENCIES
import re
import copy
import os
import random
import shutil
import multiprocessing as mp
from itertools import combinations
from itertools import repeat
from collections import Counter

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree
import numpy as np

from helper_functions import BPP_run_capture
from helper_functions import dict_to_bppcfile
from helper_functions import BPP_summary
from helper_functions import Imap_to_PopInd_Dict
from helper_functions import BPP_resume_capture
from helper_functions import path_filename

from tree_helper_functions import name_Internal_nodes
from tree_helper_functions import tree_To_Newick

from data_dicts import clprnt

from align_imap_module import autoPopParam, autoPrior


# generate a random starting newick tree
def generate_random_tree(input_node_names):
    random.shuffle(input_node_names)
    test_tree = f"({','.join(input_node_names)});"
    t = Tree(test_tree)
    t.resolve_polytomy()
    
    return tree_To_Newick(t)

# calculate the average RF distance between a set of trees
def calculate_avg_rf(tree_list):
    etetree_list = []
    for tree in tree_list:
        etetree_list.append(name_Internal_nodes(Tree(tree)))
    indexes = [i for i in range(len(tree_list))]
    tree_comparison = list(combinations(indexes, 2))
    
    rfdist = []
    for combination in tree_comparison:
        dist = etetree_list[combination[0]].compare(etetree_list[combination[1]])["rf"]
        rfdist.append(dist)

    return np.round(np.average(rfdist), decimals = 2)

def count_unique_topo(tree_list):
    ocurrences = dict(Counter(tree_list))
    return ocurrences

iteration_size = 10000    

std_cfile = {'seed': '-1',
             'outfile': 'out.txt', 
             'mcmcfile': 'mcmc.txt', 
             'speciesdelimitation': '0', 
             'speciestree': '1', 
             'usedata': '1', 
             'locusrate': '0', 
             'cleandata': '0', 
             'print': '1 0 0 0', 
             'sampfreq': '1', 
            }

# perform the iterations from the burn in up to and including the first checkpoint
def generate_tree_burinin(intree_list, imapfile, seqfile, smpl, burnin, priors, core_offset, index):
    pop_param = autoPopParam(imapfile, seqfile)
    cdict = copy.deepcopy(std_cfile)
    cdict["nsample"] = smpl
    cdict["burnin"] = burnin
    cdict["seqfile"] = f"../{seqfile}"
    cdict["imapfile"] = f"../{imapfile}"
    cdict["threads"] = f"2 {int(np.round(((index+0.5)*2), decimals = 2)+ core_offset )}"
    cdict['species&tree'] = pop_param["species&tree"]
    cdict['popsizes'] = pop_param["popsizes"]
    cdict["newick"] = intree_list[index]
    cdict['nloci'] = pop_param['nloci']
    cdict["thetaprior"] = priors['thetaprior']
    cdict["tauprior"] = priors['tauprior']
    cdict["checkpoint"] = f"{iteration_size+burnin} {iteration_size}"
    
    folder_name = f"replicate_{index}"
    os.mkdir(folder_name)
    os.chdir(folder_name)
    dict_to_bppcfile(cdict, "bpp.ctl")
    BPP_run_capture("bpp.ctl", index)
    tree = get_topology_from_MCMC()
    
    os.chdir("..")
    
    return tree

# collect the tree output of a bpp summary run
def get_topology_from_MCMC():
    full_out = BPP_summary("bpp.ctl")
    lines = full_out.split("\n")
    rowindex_tree = [i for i, s in enumerate(lines) if '(A)' in s][0]+1
    tree = re.search("\(.+\);" , lines[rowindex_tree].split("  ")[-1]).group()

    return tree

# iterate onwards from a checkpoint file and collect the next tree output
def iterate_tree_from_chk(input_folder):
    os.chdir(input_folder)
    
    ls = os.listdir()
    chk_filenames = [file for file in ls if ".chk" in str(file)]
    chk_maxval = max([int(filename.split(".")[-2]) for filename in chk_filenames])
    chk_filename = f"out.txt.{chk_maxval}.chk"
    
    BPP_resume_capture(chk_filename, input_folder.split("_")[-1])
    tree = get_topology_from_MCMC()
    os.chdir("..")
    
    return tree

# main function implementing the RF convergence testing
def test_topology(imapfile, seqfile, working_dir, repeats, smpl, burnin, core_offset = 0):
    # customized user feedback displayed in the terminal, and written to the output file
    def tree_feedback(tree_array, rf_array, samples):
        text = "\n"
        text += f"All trees after {samples} samples\n"
        for tree in tree_array[-1]: text += f"{str(tree)[1:-1]}\n"
        text += f"Unique trees and counts:\n"
        uc = count_unique_topo(tree_array[-1])
        for tree in uc: text += f"{uc[tree]} | {tree}\n"
        text += f"Average pairwise rf: {rf_array[-1]}\n"
        print(f"{clprnt.GREEN}", end = "\n")
        print(text)
        print(f"{clprnt.end}", end = "")

        return text

    # set up working directory
    parent_dir = os.getcwd()
    os.mkdir(working_dir)
    shutil.copy(src = seqfile,  dst = working_dir)
    shutil.copy(src = imapfile, dst = working_dir)
    os.chdir(working_dir)
    imapfile = path_filename(imapfile)
    seqfile = path_filename(seqfile)

    # set up output files
    with open("summary_tree_rf.csv","w") as summ_file:
        summ_file.write("average rf, samples\n")
    with open("detailed_tree_rf.txt","w") as summ_file:
        summ_file.write("")

    # set up commonly used priors
    node_names = list(Imap_to_PopInd_Dict(imapfile))
    priors = autoPrior(imapfile, seqfile)

    # generate the random starting trees that are the starting point, and initiate the RF array
    strf = 0
    while strf < 0.1:
        starting_trees = []
        # generate 2 random starting trees
        for _ in range(repeats): 
            starting_trees.append(generate_random_tree(node_names))
        strf = calculate_avg_rf(starting_trees)
    
    rf_array = [strf]
    tree_array = []
    tree_array.append(starting_trees)

    fb = tree_feedback(tree_array, rf_array, -1*burnin)
    with open("summary_tree_rf.csv","a") as summ_file:
        summ_file.write(f"{rf_array[-1]}, {-1*burnin}\n")
    with open("detailed_tree_rf.txt","a") as summ_file:
        summ_file.write(fb)


    tree_array
    # run the first iteration to start
    pool = mp.Pool(mp.cpu_count())
    treeindex = list(range(repeats))
    new_tree_list = pool.starmap(   generate_tree_burinin, 
                                    zip(repeat(tree_array[-1]), 
                                        repeat(imapfile), 
                                        repeat(seqfile), 
                                        repeat(smpl), 
                                        repeat(burnin),
                                        repeat(priors),
                                        repeat(core_offset), 
                                        treeindex))
    pool.close()

    # collect the results
    tree_array.append(new_tree_list)
    rf_array.append(calculate_avg_rf(tree_array[-1]))
    
    fb = tree_feedback(tree_array, rf_array, iteration_size)
    with open("summary_tree_rf.csv","a") as summ_file:
        summ_file.write(f"{rf_array[-1]}, {iteration_size}\n")
    with open("detailed_tree_rf.txt","a") as summ_file:
        summ_file.write(fb)

    
    # run the subsequent resume iterations
    iteration = 2
    converged = False
    folder_names = [f"replicate_{index}" for index in range(repeats)]

    while iteration*iteration_size <= smpl and converged == False:
        # run the multithreaded mode
        pool = mp.Pool(mp.cpu_count())
        treeindex = list(range(repeats))
        new_tree_list = pool.starmap(iterate_tree_from_chk, zip(folder_names))
        pool.close()

        # collect the results
        tree_array.append(new_tree_list)
        rf_array.append(calculate_avg_rf(tree_array[-1]))
        fb = tree_feedback(tree_array, rf_array, iteration_size*iteration)

        with open("summary_tree_rf.csv","a") as summ_file:
            summ_file.write(f"{rf_array[-1]}, {iteration_size*iteration}\n")
        with open("detailed_tree_rf.txt","a") as summ_file:
            summ_file.write(fb)

        # check if all of the trees have converged
        if rf_array[-1] == 0:
           converged = True 

        iteration += 1
   

    os.chdir(parent_dir)    

## EXAMPLE CALCULATIONS

test_topology   (
    imapfile="Test_Data/HLizard_2009/D_HL_imap.txt",
    seqfile="Test_Data/HLizard_2009/D_HL_align.txt", 
    working_dir="Test_Results/hliz_topo",
    repeats = 8, 
    smpl = 50000,
    burnin=20000
                )
