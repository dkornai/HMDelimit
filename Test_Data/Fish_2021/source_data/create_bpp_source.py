import os
from os import listdir
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

os.chdir("Test_Data/Fish_2021")
filenames = listdir()
filenames = [file for file in filenames if "phy" in file]
print(filenames)

# read in the individual alignments, and combine them into an MSA
full_align = []
for file in filenames:
    temp_align = AlignIO.read(file, "phylip-relaxed")
    full_align.append(temp_align)

# change the labelling of the individuals in the alignment to match with BPP formatting
for locus in full_align:
    for record in locus:
        id = record.id
        try:
            id = f"^{id.split('_')[0]}"
        except:
            id = f"^{id}"
        record.id = id
    print(locus)

# write the resulting alignment
AlignIO.write(full_align, "align_13.txt", "phylip-sequential")