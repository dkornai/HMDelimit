seed=19693
seqfile=align.txt
Imapfile=imap.txt
outfile=starting_tree_out.txt
mcmcfile=starting_tree_mcmc.txt
speciesdelimitation=0
speciestree=1
species&tree=25 aztecus_BZ aztecus_GT aztecus_HN aztecus_NI aztecus_SMX salvini_grp_1 salvini_grp_2 salvini_grp_3 salvini_grp_4 salvini_CR_a salvini_CR_b salvini_grp_5_GT_a salvini_grp_5_GT_b salvini_grp_5_HN_a salvini_grp_5_HN_b salvini_grp_5_HN_c salvini_grp_5_HN_d salvini_grp_5_HN_e salvini_grp_5_NI_a salvini_grp_6_CR_a salvini_grp_6_CR_b salvini_grp_6_CR_c salvini_CR_c salvini_PN salvini_grp_7
               2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
               (salvini_grp_2,(salvini_grp_1,((salvini_grp_3,salvini_grp_4),((salvini_grp_5_HN_a,(salvini_grp_5_HN_e,(salvini_grp_5_GT_b,(salvini_grp_5_GT_a,((salvini_grp_7,(aztecus_BZ,(aztecus_GT,(aztecus_NI,(aztecus_SMX,aztecus_HN))))),(((salvini_CR_b,salvini_CR_a),(salvini_PN,salvini_CR_c)),(salvini_grp_6_CR_a,(salvini_grp_6_CR_c,salvini_grp_6_CR_b)))))))),(salvini_grp_5_HN_b,(salvini_grp_5_NI_a,(salvini_grp_5_HN_d,salvini_grp_5_HN_c)))))));
usedata=1
nloci=50
locusrate=1 0 0 5.0 iid
cleandata=0
thetaprior=3 0.0011 e
tauprior=3 0.3914
finetune=1: 5 0.001 0.001 0.001 0.3 0.33 1.0
print=1 0 0 0
burnin=5000
sampfreq=3
nsample=15000
threads=14 2 1
