'''
THESE FUNCTIONS ARE RESPONSIBLE FOR DECIDING WHETHER THE PROPOSED CHANGES TO 
THE TOPOLOGY SHOULD BE ACCEPTED OR NOT. ADDITIONALLY, THE MODULE ALSO
WATCHES THE STOPPING CONDITIONS OF THE PIPELINE
'''
## DEPENDENCIES
# STADARD LIBRARY DEPENDENCDIES
import ast

# EXTERNAL LIBRARY DEPENDENCIES
import numpy as np

# HELPER FUNCTION DEPENDENCIES
from helper_functions import flatten
from helper_functions import readLines
from helper_functions import pretty_Table

## DATA DEPENDENCIES
from data_dicts import HM_decision_criteria_description

## TYPE HINTING DEPENDENCIES
from custom_types import BPP_out_file
from custom_types import HM_decision_parameters
from custom_types import Demographic_parameters
from custom_types import Species_name
from custom_types import Population_list
from custom_types import HM_criteria_matched


## SPECEIALIZED HELPER FUNCTIONS
# get a dict of node names:tau values and node names:theta values from the outfile
def get_Name_TauTheta_dict  (
    BPP_outfile:                    BPP_out_file
                            ) ->    tuple[dict[str, float], dict[str, float]]:

    lines = readLines(BPP_outfile)
    relevant_index = lines.index("List of nodes, taus and thetas:") # find the line where the node labels are listed
    lines = lines[relevant_index+2:]
    tau_dict = {line.split()[3]:float(line.split()[1]) for line in lines if float(line.split()[1]) != 0}
    theta_dict = {line.split()[3]:float(line.split()[2]) for line in lines}

    return tau_dict, theta_dict

# calculate the two gdi values for a given pair of species
'''
This function is used to calculate the GDI of a pair. It is core to the pipeline. 
If the theta value is not available for one of the species (this happens if only one 
sequence is assigned to that species) the gdi calculation recycles the existing
theta value for the other species, and both gdi values will be identical.
'''
def calculate_GDI   (
    Name_TauTheta_dict:     tuple[dict[str, float], dict[str, float]], 
    species_1:              Species_name, 
    species_2:              Species_name
                    ) ->    float:    

    # try to find the tau value with one or the other order of species names
    try:
        tau_key_1 = f"{species_1}{species_2}"
        tau = Name_TauTheta_dict[0][tau_key_1]
    except:
        tau_key_2 = f"{species_2}{species_1}"
        tau = Name_TauTheta_dict[0][tau_key_2]

    # try to find a theta value first for species1, but if that fails, then species 2
    try:
        theta_key_1 = f"{species_1}"
        theta = Name_TauTheta_dict[1][theta_key_1]
    except:
        theta_key_2 = f"{species_2}"
        theta = Name_TauTheta_dict[1][theta_key_2]
    
    # calculation of the gdi value
    gdi = float(np.round(1-(np.exp(-2*tau/theta)), decimals = 4))
    
    return gdi

# calculate GDI values for all elements in a list of possible mergeable/splittable pairs
def GDI_Of_Pairs(
        Name_TauTheta_dict: tuple[dict[str, float], dict[str, float]], 
        pairlist:           list[list]
                ) ->        list[list[float]]:

    output_gdis = []
    for pair in pairlist:
        output_gdis.append([calculate_GDI(Name_TauTheta_dict, pair[0], pair[1]), 
                            calculate_GDI(Name_TauTheta_dict, pair[1], pair[0])])
    
    return output_gdis

# get the most probable tau value for a given species pair
def calculate_tau   (
    Name_TauTheta_dict:     tuple[dict[str, float], dict[str, float]], 
    species_1:              Species_name, 
    species_2:              Species_name
                    ) ->    float:
    
    # try to find the tau value with one or the other order of species names
    try:
        tau_key_1 = f"{species_1}{species_2}"
        tau = Name_TauTheta_dict[0][tau_key_1]
    except:
        tau_key_2 = f"{species_2}{species_1}"
        tau = Name_TauTheta_dict[0][tau_key_2]
    
    return tau 

# calculate the # of generations since the s1 and s2 split from their common ancestor using TAU & subsitutions/site/generation
def age_Of_Pairs(
        Name_TauTheta_dict: tuple[dict[str, float], dict[str, float]], 
        pairlist:           list[list[Species_name]], 
        mutation_rate:      float
                ) ->        int:

    output_ages = []
    for pair in pairlist:
        tau = calculate_tau(Name_TauTheta_dict, pair[0], pair[1])
        age = tau/mutation_rate
        age = int(np.rint(age)) # this produces an age in whole generations
        output_ages.append(age)

    return output_ages


## FUNCTIONS IMPLEMENTING THE STEPS OF THE DECISION PROCESS

# extract the parameters (GDI, split age in generations) relevant to the merge decision from the simulation results
def get_demographic_param   (
        BPP_outfile:                BPP_out_file, 
        proposed_changes:           list[list[Species_name]], 
        hm_param:                   HM_decision_parameters,
                            ) ->    Demographic_parameters:

    name_tautheta_dict = get_Name_TauTheta_dict(BPP_outfile)
    # create empty dict to hold results
    param_dict = {str(pair):{"gdi_1": "?", "gdi_2": "?", "age": "?"} for pair in proposed_changes}
    
    # collect gdi values except if the only parameter required is "age"
    if hm_param["HM_decision"] != "age":
        gdi_values = GDI_Of_Pairs(name_tautheta_dict, proposed_changes)
        for i, pair in enumerate(proposed_changes):    
            param_dict[str(pair)]["gdi_1"] = gdi_values[i][0]
            param_dict[str(pair)]["gdi_2"] = gdi_values[i][1]
    # collect age values if the mutation rate parameter is available
    if hm_param["mutationrate"] != "?": 
        age_values = age_Of_Pairs(name_tautheta_dict, proposed_changes, hm_param["mutationrate"])
        for i, pair in enumerate(proposed_changes):    
            param_dict[str(pair)]["age"] = age_values[i]
    
    return param_dict

# based on the calculated model paramteres, and the criteria described in the master control file, decide which parameters match the required values
def criteria_matcher(
        demographic_param:  Demographic_parameters, 
        hm_param:           HM_decision_parameters
                    ) ->    HM_criteria_matched:

    # create empty result dict
    match_dict = {str(pair):{"gdi_1": "?", "gdi_2": "?", "age": "?"} for pair in demographic_param}

    age_thresh = hm_param["generations"]
    GDI_thresh = hm_param["GDI_thresh"]

    # iterate through all the possible pairs and decide which criteria they match
    for pair in demographic_param.keys():
        # seperate the variables out for shorter conditional checks
        gdi_1 = demographic_param[pair]["gdi_1"]
        gdi_2 = demographic_param[pair]["gdi_2"]
        age   = demographic_param[pair]["age"]
    
        # if gdi values are available, decide if the gdi value is above or below the threshold 
        if gdi_1 != "?" and gdi_2 != "?":
            if hm_param["mode"] == "merge":
                if gdi_1 <= GDI_thresh:
                    match_dict[str(pair)]["gdi_1"] = True
                else:
                    match_dict[str(pair)]["gdi_1"] = False
                if gdi_2 <= GDI_thresh:
                    match_dict[str(pair)]["gdi_2"] = True
                else:
                    match_dict[str(pair)]["gdi_2"] = False
            elif hm_param["mode"] == "split":
                if gdi_1 >= GDI_thresh:
                    match_dict[str(pair)]["gdi_1"] = True
                else:
                    match_dict[str(pair)]["gdi_1"] = False
                if gdi_2 >= GDI_thresh:
                    match_dict[str(pair)]["gdi_2"] = True
                else:
                    match_dict[str(pair)]["gdi_2"] = False
        
        # if an age value is available, decide if the age of the splits is above or below the threshold
        if age != "?":
            if hm_param["mode"] == "merge":
                if age <= age_thresh:
                    match_dict[str(pair)]["age"] = True
                else:
                    match_dict[str(pair)]["age"] = False
            elif hm_param["mode"] == "split":
                if age >= age_thresh:
                    match_dict[str(pair)]["age"] = True
                else:
                    match_dict[str(pair)]["age"] = False

    return match_dict

# return the list of proposal pairs that match the correct number and type of parameters to be accepted
def make_decision   (
        input_match_dict:   HM_criteria_matched, 
        hm_param:           HM_decision_parameters
                    ) ->    list[list[Species_name]]:

    # set up empty container to store results
    accepted = []
    # evaulate incoming statements depending on the decision criteria
    for pair in input_match_dict:
        pair_as_list = ast.literal_eval(pair)

        # seperate the variables out for shorter conditional checks
        gdi_1 = input_match_dict[pair]["gdi_1"]
        gdi_2 = input_match_dict[pair]["gdi_2"]
        age   = input_match_dict[pair]["age"]
        results = [gdi_1, gdi_2, age]

        # all proposals are accepted regardless of results, this is useful for exploring the GDI values 
        if hm_param["HM_decision"]   == "none":         
            accepted.append(pair_as_list)
        # a minimum of one of the parameters within the threshold to accept a proposal
        elif hm_param["HM_decision"] == "any":          
            if results.count(True) > 0:
                accepted.append(pair_as_list)
        # two of the parameters must be within their threshold to accept a proposal
        elif hm_param["HM_decision"] == "any_two":        
            if results.count(True) > 1:
                accepted.append(pair_as_list)
        # only one gdi is enough to accept a proposal
        elif hm_param["HM_decision"] == "one_gdi":      
            if gdi_1 == True or gdi_2 == True:
                accepted.append(pair_as_list)
        # both gdis need to be within their threshold to accept a proposal
        elif hm_param["HM_decision"] == "both_gdis":    
            if gdi_1 == True and gdi_2 == True:
                accepted.append(pair_as_list)
        # only age is enough to accept a proposal
        elif hm_param["HM_decision"] == "age":                 
            if age == True:
                accepted.append(pair_as_list)
        # age and one gdi needs to be within the threshold to accept a proposal
        elif hm_param["HM_decision"] == "one_gdi_&_age":
            if (gdi_1 == True or gdi_2 == True) and age == True:
                accepted.append(pair_as_list)
        # all calculated parameters must be within their threshold to accept a proposal
        elif hm_param["HM_decision"] == "all": 
            if results.count(True) == (len(results) - results.count("?")):
                accepted.append(pair_as_list)           

    return accepted

# prints feedback to the user about the decision criteria
def decisionUserFeedback(
        proposed_changes:       list[list[Species_name]],
        demographic_parameters: Demographic_parameters, 
        hm_param:               HM_decision_parameters, 
        matched_dict:           HM_criteria_matched,   
        decision:               list[list[Species_name]]
                        ):

    # collect if the program is working in merge or split mode
    mode = hm_param['mode'].upper()
    gdi_verb = "<"
    if mode == "SPLIT":
        gdi_verb = ">"
    
    # feedback about the threshold used in the decisions
    print(f"\n1) The {mode} decision thresholds are:")
    if hm_param['HM_decision'] != "age":
        print(f"\n  GDI values must be: {gdi_verb} {hm_param['GDI_thresh']} to be considered sufficient")
    if hm_param['mutationrate'] != "?":
        print(f"\n  The populations should share a common ancestor: {gdi_verb} {hm_param['generations']} generations ago")

    # feedback about the criteria used in the decisions
    print(f"\n2) The decision criteria is: {hm_param['HM_decision'].upper()}")
    print(f"\n  {mode.lower()} proposals are accepted {HM_decision_criteria_description[hm_param['HM_decision']]}")
    
    # gather the decision demographic_parameters and results into list of lists
        
        # names of the population pairs that are proposed to change
    name_list = []
    for proposal in proposed_changes:
        name_list.append(str(proposal))
    results_table = [name_list]
    results_colnames = [f"POPULATIONS TO {mode}"]
        
        # if GDIs are calculated, collect GDI values, and GDI acceptances
    if hm_param["HM_decision"] != "age":
        gdi_1_list = []
        gdi_2_list = []
        gdi_1_acc_list = []
        gdi_2_acc_list = []

        for proposal in proposed_changes:
            gdi_1_list.append(str(demographic_parameters[str(proposal)]["gdi_1"]))
            gdi_1_acc_list.append(str(matched_dict[str(proposal)]["gdi_1"]))
            gdi_2_list.append(str(demographic_parameters[str(proposal)]["gdi_2"]))
            gdi_2_acc_list.append(str(matched_dict[str(proposal)]["gdi_2"]))

        results_table.extend([gdi_1_list, gdi_1_acc_list , gdi_2_list, gdi_2_acc_list])
        results_colnames.extend(["GDI 1", "SUFFICIENT", "GDI 2", "SUFFICIENT"])

        # if age is calculated, collect age values and age acceptances
    if hm_param['mutationrate'] != "?":
        age_list = []
        age_acc_list = []

        for proposal in proposed_changes:
            age_list.append(str(demographic_parameters[str(proposal)]["age"]))
            age_acc_list.append(str(matched_dict[str(proposal)]["age"]))
        
        results_table.extend([age_list, age_acc_list])
        results_colnames.extend(["# OF GENERATIONS", "SUFFICIENT"])

        # collect if the proposasl was accepted or not
    accepted_list = []
    for proposal in proposed_changes:   
        accepted = False
        if proposal in decision:
            accepted = True
        accepted_list.append(str(accepted))
    results_table.extend([accepted_list])
    results_colnames.extend([f"{mode} ACCEPTED"])
    
    # print a table detailing all results to the user
    print("\nTABLE OF RESULTS:\n")
    pretty_Table(results_table, results_colnames)
    print()

    # print the list of accepted proposals
    if len(decision) == 0:
        print(f"NO PROPOSALS TO {mode} WERE ACCEPTED")
    else:
        print(f"THE FOLLOWING PROPOSALS TO {mode} WERE ACCEPTED:")
        for i, pair in enumerate(decision):
            print(f"{i+1}) {str(pair)[1:-2]}")

    

# implements the decisions made by previous modules by changing the list of accepted populations
def implement_decision  (
        previous_pops:          Population_list, 
        decision:               list[list[Species_name]],
        hm_param:               HM_decision_parameters, 
                        ) ->    Population_list:

    flattened_decision_list:Population_list = flatten(decision)

    # in merge mode, the population that meet the criteria are removed (merged into their ancestor, which remains)
    if hm_param["mode"] == "merge":
        new_accepted_pops = [pop for pop in previous_pops if pop not in flattened_decision_list]
    
    # in split mode, the population that meet the criteria are added (split from their ancestor, which also remains)
    elif hm_param["mode"] == "split":
        new_accepted_pops = previous_pops + flattened_decision_list
    
    print("\nACCORDINGLY, THE LIST OF POPULATIONS AND ANCESTRAL POPULATIONS THAT ARE CURRENTLY ACCEPTED AS SPECIES IS:\n")
    for population in new_accepted_pops:
        print(f"\t{population}")
    print()

    return new_accepted_pops


# check if the Hierarchical merge process has reached any of the end conditions, and print feedback to the user
def stop_check  (
        hm_param:           HM_decision_parameters, 
        decision:           list[list[Species_name]], 
        new_accepted_pops:  Population_list, 
        halt_pop_number:    int
                ) ->        bool:

    # set feedback text according to mode
    mode = hm_param['mode'].upper()
    if mode == "MERGE":
        end_topo = "CONSISTS OF JUST THE ROOT NODE OF THE GUDIE TREE"
    else:
        end_topo = "HAS REACHED THE GUIDE TREE"

    #    A) no new proposals were accepted
    if len(decision) == 0:
        to_iterate = False
        
        print(f"AS ALL {mode} PROPOSALS WERE REJECTED, NO FURTHER PROPOSALS CAN BE MADE!")
    
    #    B) the delimitation has reached either the root, or the guide tree, and no further moves can be made
    elif len(new_accepted_pops) == halt_pop_number:
        to_iterate = False
        
        print(f"AS THE ACCEPTED TOPOLOGY {end_topo}, NO FURTHER {mode} PROPOSALS CAN BE MADE!")

    # if none of the end conditions were reached, continue the program
    else:
        to_iterate = True

    return to_iterate


## FINAL WRAPPER FUNCTION IMPLEMENTING THE COMPLETE DECISION PROCESS

# wrapper function that implements the complete decision procedure
def decisionModule  (
        hm_param:           HM_decision_parameters, 
        BPP_outfile:        BPP_out_file, 
        proposed_changes:   list[list[Species_name]], 
        accepted_pops:      Population_list, 
        halt_pop_number:    int
                    ) ->    tuple[Population_list, bool]:

    print("\nMAKING DECISIONS BASED ON BPP MODEL RESULTS")

    # extract the parameter values relevant to the decision
    demog_param = get_demographic_param(BPP_outfile, proposed_changes, hm_param)
    
    # check if the parameter values are within the thresholds required to accept a decision
    match_dict = criteria_matcher(demog_param, hm_param)
    
    # get the final list of population pairs that match the necessary criteria to accept
    decision = make_decision(match_dict, hm_param)

    # print feedback to the user about the decision process
    decisionUserFeedback(proposed_changes, demog_param, hm_param, match_dict, decision)

    # implement the changes to the list of accepted popuations
    new_accepted_pops = implement_decision(accepted_pops, decision, hm_param)

    # keep track of whether the program has finished 
    to_iterate = stop_check(hm_param, decision, new_accepted_pops, halt_pop_number)

    return new_accepted_pops, to_iterate