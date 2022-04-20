          seed =  -1

       seqfile = ../lizard.txt
      Imapfile = ../lizard.Imap.txt
       outfile = out
      mcmcfile = mcmc.out

 speciesdelimitation = 0 * fixed species tree
 speciesdelimitation = 1 0 5    * speciesdelimitation algorithm0 and finetune(e)
*  speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)
   speciestree = 1

  species&tree = 5  tri  cow  con  und  woo
                     4    3    4    5    1
                 (((tri, cow), con), (und, woo));

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 1 29    * number of data sets in seqfile

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 1000    # gamma(a, b) for theta
      tauprior = 2 1000 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

       finetune = 1: 5 0.0005 0.002  0.0005 0.5 0.2 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 0
        burnin = 5000
      sampfreq = 10
       nsample = 100000

*** Note: Make your window wider (140 columns) before running the program.
