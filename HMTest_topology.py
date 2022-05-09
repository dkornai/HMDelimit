
from unittest import result
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree

from helper_functions import BPP_run, Imap_to_PopInd_Dict
from helper_functions import extract_Speciestree
from helper_functions import dict_to_bppcfile
import copy
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import random
from itertools import combinations
from itertools import repeat
from tree_helper_functions import name_Internal_nodes
from tree_helper_functions import tree_To_Newick
from align_imap_module import autoPopParam, autoPrior
import shutil
from helper_functions import path_filename

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


def tree_iterate(intree_list, imapfile, seqfile, smpl_interval, iteration, priors, index):
    pop_param = autoPopParam(imapfile, seqfile)
    cdict = copy.deepcopy(std_cfile)
    cdict["nsample"] = smpl_interval
    cdict["burnin"] = 0
    cdict["seqfile"] = f"../{seqfile}"
    cdict["imapfile"] = f"../{imapfile}"
    cdict["threads"] = f"1 {int(index+1)}"
    cdict['species&tree'] = pop_param["species&tree"]
    cdict['popsizes'] = pop_param["popsizes"]
    cdict["newick"] = intree_list[index]
    cdict['nloci'] = pop_param['nloci']
    cdict["thetaprior"] = priors['thetaprior']
    cdict["tauprior"] = priors['tauprior']
    folder_name = f"smpl_{smpl_interval*iteration}_{index}.txt"
    os.mkdir(folder_name)
    os.chdir(folder_name)
    dict_to_bppcfile(cdict, "bpp.ctl")
    print(os.getcwd())
    BPP_run("bpp.ctl")
    output_tree = extract_Speciestree("bpp.ctl")
    os.chdir("..")
    
    return output_tree



def test_topology(imapfile, seqfile, working_dir, repeats, smpl_interval, max_sample):
    parent_dir = os.getcwd()
    os.mkdir(working_dir)
    shutil.copy(src = seqfile,  dst = working_dir)
    shutil.copy(src = imapfile, dst = working_dir)
    os.chdir(working_dir)
    imapfile = path_filename(imapfile)
    seqfile = path_filename(seqfile)

    r_smpl = []
    r_sites = []
    r_rf_s = []
    r_rf_e = []
    r_rf_r = []

    node_names = list(Imap_to_PopInd_Dict(imapfile))
    priors = autoPrior(imapfile, seqfile)

    # generate the random starting trees that are the starting point
    strf = 0
    while strf < 0.1:
        starting_trees = []
        # generate 2 random starting trees
        for _ in range(repeats): 
            starting_trees.append(generate_random_tree(node_names))
        strf = calculate_avg_rf(starting_trees)
    
    print("@@@@@ STARTING TREES @@@@@@@")
    for t in starting_trees: print(t)

    # begin the tree array with the starting trees
    tree_array = []
    tree_array.append(starting_trees)
    rf_array = [strf]

    iteration = 1
    while iteration*smpl_interval <= max_sample:
        
        # run the data generation in multithreaded mode
        pool = mp.Pool(mp.cpu_count())
        treeindex = list(range(repeats))

        new_tree_list = pool.starmap(tree_iterate, zip(repeat(tree_array[-1]), 
                                                       repeat(imapfile), 
                                                       repeat(seqfile), 
                                                       repeat(smpl_interval), 
                                                       repeat(iteration),
                                                       repeat(priors), treeindex))

        pool.close()    
        tree_array.append(new_tree_list)
        rf_array.append(calculate_avg_rf(new_tree_list))
        print(f"@@@@@ TREES AFTER ITERATION {iteration} @@@@@@@")
        for t in new_tree_list: print(t)
        iteration += 1
        print("RF:", rf_array[-1])
            
        # start_rf = calculate_avg_rf(starting_trees)
        # end_rf = calculate_avg_rf(end_trees)
        # rf_ratio = np.round(end_rf/start_rf, decimals = 2)

        # r_smpl.append(sample)
        # r_rf_s.append(start_rf)
        # r_rf_e.append(end_rf)
        # r_rf_r.append(rf_ratio)

   
    results = pd.DataFrame({"rf":rf_array,})
    results.to_csv(path_or_buf=f"results.csv", index=False)
    os.chdir(parent_dir)    
    return rf_array, tree_array

rf_vals, trees = test_topology("Test_Data/Snakes_2019/imap.txt",
              "Test_Data/Snakes_2019/align_11.txt", 
              "Test_Results/snk",
              18, 10000, 250000)

rf_vals, trees = test_topology("Test_Data/Fish_2021/imap.txt",
              "Test_Data/Fish_2021/align_13.txt", 
              "Test_Results/fsh",
              18, 10000, 200000)

rf_vals, trees = test_topology("Test_Data/HLizard_2009/D_HL_imap.txt",
              "Test_Data/HLizard_2009/D_HL_align.txt", 
              "Test_Results/hliz",
              18, 10000, 250000)

rf_vals, trees = test_topology("Test_Data/Sarracenia_2013/imap.txt",
              "Test_Data/Sarracenia_2013/align_a.txt", 
              "Test_Results/sarr",
              18, 10000, 250000)

rf_vals, trees = test_topology("Test_Data/TMS_2019/D_TMS_imap.txt",
              "Test_Data/TMS_2019/D_TMS_align.txt", 
              "Test_Results/tms",
              18, 10000, 250000)
