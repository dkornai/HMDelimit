seed=20365
seqfile=align.txt
Imapfile=Imap_UniqueID.txt
outfile=starting_delimitation_out.txt
mcmcfile=starting_delimitation_mcmc.txt
speciesdelimitation=1 1 2 1
speciestree=1
species&tree=25 PN_000 PN_001 PN_002 PN_003 PN_004 PN_005 PN_006 PN_007 PN_008 PN_009 PN_010 PN_011 PN_012 PN_013 PN_014 PN_015 PN_016 PN_017 PN_018 PN_019 PN_020 PN_021 PN_022 PN_023 PN_024  
               2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
               (((((((((PN_000,PN_001),PN_002),PN_003),PN_004),PN_024),((((PN_009,(PN_010,PN_022)),PN_023),(PN_019,(PN_020,PN_021))),((PN_011,PN_012),(((PN_013,PN_014),((PN_015,PN_018),PN_016)),PN_017)))),PN_008),PN_007),(PN_005,PN_006));
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
