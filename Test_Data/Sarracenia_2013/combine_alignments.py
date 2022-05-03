import copy
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from helper_functions import alignfile_to_MSA

def report_informative(input_MSA:MultipleSeqAlignment):
    length = input_MSA.get_alignment_length()
    print("length is:",length)
    
    informative = 0
    for i in range(length):
        if len(set(str(input_MSA[:, i]))) > 1:
            informative += 1

    print("the number of informative sites is:", informative)
    return informative


align_a = alignfile_to_MSA("align_a.txt")
# measure the number of informative sites
info_1 = 0    
for locus in align_a:
    info_1 += report_informative(locus)
print(info_1)
align_b = alignfile_to_MSA("align_b.txt")
# measure the number of informative sites
info_1 = 0    
for locus in align_b:
    info_1 += report_informative(locus)
print(info_1)

# create combined alignment
combined_align = []

for i in range(len(align_a)):
    temp_locus = copy.deepcopy(align_a[i])
    temp_locus.extend(copy.deepcopy(align_b[i]))
    combined_align.append(temp_locus)

# measure the number of informative sites
info_2 = 0    
for locus in combined_align:
    info_2 += report_informative(locus)
print(info_2)

AlignIO.write(combined_align, "align_combined.txt", "phylip-sequential")