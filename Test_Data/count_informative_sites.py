
from helper_functions import alignfile_to_MSA

aln = alignfile_to_MSA("Test_Data/CLizard_2017/align.txt")

def report_informative(input_MSA):
    length = input_MSA.get_alignment_length()
    print("length is:",length)
    
    informative = 0
    for i in range(length):
        if len(set(str(input_MSA[:, i]))) > 1:
            informative += 1

    print("the number of informative sites is:", informative)
    return informative

# measure the number of informative sites
info = 0    
for locus in aln:
    info += report_informative(locus)
print(info)