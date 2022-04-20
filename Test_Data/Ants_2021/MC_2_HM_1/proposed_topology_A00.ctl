seed=35156
seqfile=align.txt
Imapfile=proposed_imap.txt
outfile=proposal_out.txt
mcmcfile=proposal_mcmc.txt
speciesdelimitation=0
speciestree=0
species&tree=20 aztecus_BZaztecus_GT aztecus_HN aztecus_NI aztecus_SMX salvini_grp_1 salvini_grp_2 salvini_grp_3 salvini_grp_4 salvini_CR_a salvini_CR_b salvini_grp_5_GT_asalvini_grp_5_GT_b salvini_grp_5_HN_a salvini_grp_5_HN_b salvini_grp_5_HN_csalvini_grp_5_NI_a salvini_grp_5_HN_d salvini_grp_5_HN_e salvini_grp_6_CR_asalvini_grp_6_CR_bsalvini_grp_6_CR_c salvini_CR_c salvini_PN salvini_grp_7
               4  2  2  2  2  2  2  2  2  2  4  2  2  4  2  2  6  2  2  2
               (((((((aztecus_BZaztecus_GT,(aztecus_HN,aztecus_NI)),aztecus_SMX),salvini_grp_7),((((salvini_CR_a,(salvini_CR_b,salvini_CR_c)),salvini_PN),salvini_grp_6_CR_asalvini_grp_6_CR_bsalvini_grp_6_CR_c),(salvini_grp_5_GT_asalvini_grp_5_GT_b,(((salvini_grp_5_HN_a,salvini_grp_5_HN_b),(salvini_grp_5_HN_csalvini_grp_5_NI_a,salvini_grp_5_HN_d)),salvini_grp_5_HN_e)))),salvini_grp_4),salvini_grp_3),(salvini_grp_1,salvini_grp_2));
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
