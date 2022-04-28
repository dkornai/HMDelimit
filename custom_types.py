'''
THESE CUSTOM TYPES ARE USED AS TYPE HINTS TO MAKE 
CODE MORE UNDERSTANDABLE
'''
from typing import NewType

# generic filepath type
file_path = NewType("file_path", str)

# generic type for text files that have been converted into a list of rows
Text_rows_list = NewType("Text_rows_list", list)

# custom types for BPP control files, and control dicts.
BPP_control_file = NewType("BPP_control_file", file_path)
BPP_control_dict = NewType("BPP_control_dict", dict)
BPP_control_dict_component = NewType("BPP_control_dict_component", dict)
BPP_out_file = NewType("BPP_out_file", file_path)

# custom types for Master control files, and associated control dicts.
Master_control_file = NewType("Master_control_file", file_path)
Master_control_dict = NewType("Master_control_dict", dict)

# custom types associated with Imap files
Imap_file = NewType("Imap_file", file_path)
Imap_list = NewType("Imap_list", list)

# custom types for a newick tree
Tree_newick = NewType("Tree_newick", str)

# custom type for a phylip alignment
Phylip_MSA_file = NewType("Phylip_MSA_file", file_path)

# custom type for a lists of population names
Population_list = NewType("Population_list", list)

# custom type for the alias dicts used when match a population name with some other version of itself, eg in the "uniqueID" encoding and decoding process
Alias_dict = NewType("Alias_dict", dict)

# custom type for BPP modes (A00, A01, A11)
BPP_mode = NewType("BPP_mode", str)

# custom type for the mode of the HM pipeline (merge or split)
HM_mode = NewType("HM_mode", str)

# custom type for species names
Species_name = NewType("Species_name", str)

# custom types for the decision module
HM_decision_parameters = NewType("HM_decision_parameters", dict)
MSC_parameters = NewType("MSC_parameters", dict)
HM_criteria_matched = NewType("HM_criteria_matched", dict)