seed = -1

seqfile = coding_pps_50_best_partitions_bpp_seqs.txt
Imapfile = coding_pps_50_best_partitions_bpp_map.txt
outfile = out.txt
mcmcfile = mcmc.txt

speciesdelimitation = 1 1 2 1
*speciestree = 1


species&tree = 25  aztecus_BZ aztecus_BZ aztecus_HN aztecus_NI aztecus_SMX salvini_grp_1 salvini_grp_2 salvini_grp_3 salvini_grp_4 salvini_CR_a salvini_CR_b salvini_grp_5_GT_a salvini_grp_5_GT_b salvini_grp_5_HN_a salvini_grp_5_HN_b salvini_grp_5_HN_c salvini_grp_5_HN_d salvini_grp_5_HN_e salvini_grp_5_NI_a salvini_grp_6_CR_a salvini_grp_6_CR_b salvini_grp_6_CR_c salvini_CR_c salvini_PN salvini_grp_7
                   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
                   ((salvini_grp_1,salvini_grp_2),(salvini_grp_3,(salvini_grp_4,((salvini_grp_7,((aztecus_HN,aztecus_NI),(aztecus_GT,(aztecus_SMX,aztecus_BZ)))),(((salvini_grp_6_CR_c,(salvini_grp_6_CR_a,salvini_grp_6_CR_b)),(salvini_PN,(salvini_CR_a,(salvini_CR_c,salvini_CR_b)))),((salvini_grp_5_GT_b,salvini_grp_5_GT_a),(salvini_grp_5_HN_e,((salvini_grp_5_HN_c,salvini_grp_5_NI_a),(salvini_grp_5_HN_d,(salvini_grp_5_HN_a,salvini_grp_5_HN_b))))))))));
                   diploid =  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

usedata = 1
nloci = 50
locusrate = 1 5
cleandata = 0
thetaprior = 3 0.06
tauprior = 3 0.01
finetune = 1: 5 0.001 0.001 0.001 0.3 0.33 1.0
print = 1 0 0 0
burnin = 50000
sampfreq = 5
nsample = 200000

mcfile

sequencefile =

locusraaate    =  

foobar