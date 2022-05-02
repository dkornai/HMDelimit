          seed = -1

       seqfile = ../RotariaSub.txt
      Imapfile = ../Rotaria.Imap.txt
       outfile = out
      mcmcfile = mcmc.out

 speciesdelimitation = 0 * fixed species tree
 speciesdelimitation = 1 0 2    * speciesdelimitation algorithm0 and finetune(e)
 speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)
        speciestree = 1

  species&tree = 4  T  N  S  M
                    8  8  12 8
                 (((T, N), S), M);

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 2    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?


    thetaprior = 2 40    # gamma(a, b) for theta
      tauprior = 2 200 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

      heredity = 0 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)

     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)

 sequenceerror = 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

      finetune = 1: 0.3  0.0005  0.05  0.001  0.06  0.25  1.0 # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 0
        burnin = 10000
      sampfreq = 2
       nsample = 100000

*** Note: Make your window wider (120 columns) before running the program.

