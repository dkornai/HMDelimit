
# from checkinput_module import check_Finetune, check_Numeric, check_Tauprior, check_Thetaprior, imap_To_Dict, pretty, proposal
# from data_dicts import HM_decision_criteria

# print(check_Finetune("1: .01 .0001 .005 .0005 .2 .01 .01 .01"),
#       check_Finetune("1: 0.01 0.004 0.0003 0.0001 .02"),
#       check_Finetune("1 0.01 0.004 0.0003 0.0001 .02"),
#       check_Finetune("2: 0.01 0.004 0.0003 0.0001 .02"),
#       check_Finetune("2: 0.01 3 0.0003 0.0001 .02"),
#       )

# print(check_Tauprior("3 0.004"),
#       check_Tauprior("3 0.004 4"),
#       check_Tauprior("3 3"))

# print(check_Thetaprior("3 0.004 e"),
#       check_Thetaprior("3 0.004 "),
#       check_Thetaprior("3 0.004 e"),)

# print(check_Numeric("3", "f"), 
#       check_Numeric("3", "i"),
#       check_Numeric("3", "f", min = 2, max = 50),
#       check_Numeric("0.3", "i"),
#       check_Numeric("0.3", "i", min= 0, max= 1),
#       check_Numeric("0.3", "i", min= 0, max= 0.2),
#       check_Numeric("?", "i", min= 0, max= 0.2),
#       check_Numeric("asdasd", "i", min= 0, max= 0.2),
# )

# # test of criteriaMatch
# test_dict_for_criteria =   {"+G1+G2+A":{"gdi_1": True,  "gdi_2": True,  "age": True},
#                             "+G1+G2-A":{"gdi_1": True,  "gdi_2": True,  "age": False},
#                             "+G1-G2+A":{"gdi_1": True,  "gdi_2": False, "age": True},
#                             "-G1+G2+A":{"gdi_1": False, "gdi_2": True,  "age": True},
#                             "+G1-G2-A":{"gdi_1": True,  "gdi_2": False, "age": False},
#                             "-G1+G2-A":{"gdi_1": False, "gdi_2": True,  "age": False},
#                             "-G1-G2+A":{"gdi_1": False, "gdi_2": False, "age": True},
#                             "-G1-G2-A":{"gdi_1": False, "gdi_2": False, "age": False},
#                             "+G1+G2?A":{"gdi_1": True,  "gdi_2": True,  "age": "?"},
#                             "+G1-G2?A":{"gdi_1": True,  "gdi_2": False, "age": "?"},
#                             "-G1-G2?A":{"gdi_1": False, "gdi_2":False,  "age": "?"},
#                             "?G1?G2+A":{"gdi_1": "?",   "gdi_2": "?",   "age": True},
#                             "?G1?G2-A":{"gdi_1": "?",   "gdi_2": "?",   "age": False},
#                }
# for criteria in HM_decision_criteria:
#     print(criteria, criteriaMatch(test_dict_for_criteria, hm_param = {"HM_decision":criteria}))

# # pretty_Table testing
# from helper_functions import pretty_Table
# column_1 = ["1", "23", "456", "0", "?"]
# column_2 = ["asd", "asdasd", "asdaa", "asd", "a"]
# column_3 = ["True", "False", "True", "False", "?"]
# table = [column_1, column_2, column_3]
# colnames = ["NUMBER", "LETTER", "BOOL"]
# pretty_Table(table, colnames)

# from proposal_module import HMproposal
# from helper_functions import Imap_to_IndPop_Dict
# # generate a proposal
# proposed = HMproposal(guide_tree_newick = "(((CBC, (NBC, SCA)), NCA), SBC);",
#                       base_indpop_dict    = Imap_to_IndPop_Dict("D_HL_imap.txt"),
#                       current_pops_list = ['CBCNBCSCANCASBC', 'CBCNBCSCANCA', 'SBC'],#, 'CBCNBCSCA', 'NCA'],
#                       mode              = "merge")
# print(proposed)

# # decisionModule testing
# from decision_module import get_HM_parameters
# hm_par = get_HM_parameters("Pcontrol.txt")

# from decision_module import decisionModule
# decisionModule("proposal_mcmc.txt", [['CBCNBCSCA', 'NCA']], "Pcontrol.txt")

# # tree maker testing
# from align_imap_module import generateStartingTree

# t = generateStartingTree("D_HL_imap.txt", "D_HL_align.txt")
# print(t)

# check the main HM iteration
# from helper_functions import Imap_to_IndPop_Dict
# from stage_modules import HMIteration
# HMIteration("Pcontrol.txt", 
#             "(((CBC, (NBC, SCA)), NCA), SBC);", 
#             Imap_to_IndPop_Dict("D_HL_imap.txt"), 
#             ['CBCNBCSCANCASBC', 'CBCNBCSCANCA', 'SBC', 'CBCNBCSCA', 'NCA'], 
#             1)

# HMIteration("Pcontrol2.txt", 
#             "(((CBC, (NBC, SCA)), NCA), SBC);", 
#             Imap_to_IndPop_Dict("D_HL_imap.txt"), 
#             ['CBCNBCSCANCASBC', 'CBCNBCSCANCA', 'SBC', 'CBCNBCSCA', 'NCA'], 
#             2)

# HMIteration("Pcontrol3.txt", 
#             "(((CBC, (NBC, SCA)), NCA), SBC);", 
#             Imap_to_IndPop_Dict("D_HL_imap.txt"), 
#             ['CBCNBCSCANCASBC', 'CBCNBCSCANCA', 'SBC', 'CBCNBCSCA', 'NCA'], 
#             3)

# HMIteration("Pcontrol4.txt", 
#             "(((CBC, (NBC, SCA)), NCA), SBC);", 
#             Imap_to_IndPop_Dict("D_HL_imap.txt"), 
#             ['CBCNBCSCANCASBC', 'CBCNBCSCANCA', 'SBC', 'CBCNBCSCA', 'NCA'], 
#            4)

# # check input testing
# from checkinput_module import check_Master_Control
# check_Master_Control("Pcontrol.txt")
# check_Master_Control("MC_badparam.txt")

# # checking of imap and sequence compatibility checker
# from align_imap_module import check_Imap_Seq_compat
# check_Imap_Seq_compat("D_HL_imap.txt", "D_HL_align.txt")

# check_Imap_Seq_compat("D_HL_imap.txt", "D_TMS_align.txt")

# # checking of imap and tree compatibility checker
# from align_imap_module import check_Imap_Tree_compat
# check_Imap_Tree_compat("D_HL_imap.txt", "(((CBC, (NBC, SCA)), NCA), SBC);")
# check_Imap_Tree_compat("D_HL_imap.txt", "(((A, B), (C, D)), X);")

# from checkinput_module import check_Alignment
# print(check_Alignment("D_HL_align.txt"))
# print(check_Alignment("D_TMS_align.txt"))

# # checking of the autoprior module
# from align_imap_module import autoPrior
# print("HL")
# print(autoPrior("D_HL_imap.txt", "D_HL_align.txt"))
# print("TMS")
# print(autoPrior("D_TMS_imap.txt", "D_TMS_align.txt"))
# print("ROT")
# print(autoPrior("D_ROT_imap.txt", "D_ROT_align.txt"))
# print("L10")
# print(autoPrior("D_L10_imap.txt", "D_L10_align.txt"))

# # checking of the starting tree module
# from align_imap_module import autoStartingTree
# print("HL")
# print(autoStartingTree("D_HL_imap.txt", "D_HL_align.txt"))
# print("TMS")
# print(autoStartingTree("D_TMS_imap.txt", "D_TMS_align.txt"))
# print("ROT")
# print(autoStartingTree("D_ROT_imap.txt", "D_ROT_align.txt"))
# print("L10")
# print(autoStartingTree("D_L10_imap.txt", "D_L10_align.txt"))

# # checking of the Imap-MSA and Imap-Tree compatibility checkers
# from align_imap_module import check_Imap_Tree_compat
# from align_imap_module import check_Imap_Seq_compat
# check_Imap_Seq_compat("D_ROT_imap.txt", "D_ROT_align.txt")
# check_Imap_Seq_compat("D_ROT_imap.txt", "D_L10_align.txt")
# check_Imap_Tree_compat("D_ROT_imap.txt", "(cow,(woo,(tri,(und,con))));")
# check_Imap_Tree_compat("D_L10_imap.txt", "(cow,(woo,(tri,(und,con))));")

# from align_imap_module import check_GuideTree_Imap_MSA_compat
# check_GuideTree_Imap_MSA_compat("(((D,C),(B,A)),X);", "D_TMS_imap.txt", "D_TMS_align_mod.txt")
# check_GuideTree_Imap_MSA_compat("(((D,C),(B,A)),X);", "D_TMS_imap.txt", "D_TMS_align.txt")
# check_GuideTree_Imap_MSA_compat("((S,(N,T)),M);", "D_ROT_imap.txt", "D_ROT_align.txt")

# from stage_modules import StartingDelimitation
# tree, pops = StartingDelimitation("Pcontrol.txt")
# print(tree)
# print(pops)