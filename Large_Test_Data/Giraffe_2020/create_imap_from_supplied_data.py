from helper_functions import alignfile_to_MSA
from helper_functions import readLines
from helper_functions import list_To_Imap

gir = alignfile_to_MSA("BPP_Giraffa66.txt")

ind_ids = []
for locus in gir:
    for seq in locus:
        ind_ids.append(seq.id)

ind_ids = list(set(ind_ids))
ind_ids = [indid.split("^")[1] for indid in ind_ids]
print(ind_ids)
print(len(ind_ids), "\n")


species_map = readLines("Input_STRUCTURE_noassignation.txt")
species_map = species_map[1:]
species_map = [line.split()[0] for line in species_map]
species_map = list(set(species_map))
print(species_map, "\n")

imap = {}
for indid in ind_ids:
    for speciesid in species_map:
        if indid in speciesid:
            imap[indid] = speciesid.split("_")[0]
imap = {k: v for k, v in sorted(imap.items(), key=lambda item: item[1])}

print(imap, "\n")

imaplist = [list(imap.keys()), list(imap.values())]

print(imaplist, "\n")

list_To_Imap(imaplist, "imap.txt")