          seed =  -1

       seqfile = D_L10_align.txt
      Imapfile = D_L10_imap.txt

  species&tree = 5  tri  cow  con  und  woo
                     4    3    4    5    1
                 (((tri, cow), con), (und, woo));

       usedata = 1 
         nloci = 29 

     cleandata = 1    

    thetaprior = 3 0.012 e   
      tauprior = 3 0.1 

       finetune = 1: 5 0.0005 0.002  0.0005 0.5 0.2 1.0 

        burnin = 10000
      sampfreq = 5
       nsample = 10000    

	threads = 8
	
*** Note: Make your window wider (140 columns) before running the program.
