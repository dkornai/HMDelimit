# HMDelimit
This repo holds the python source files for a species delimitation pipeline based on the work of Leache et al 2019.
# Methods:
The pipeline performs multilocus species delimitation by combining the use of BPP, and a metric of speciation called the Genomic Divergence Index(GDI) 
### A simple example:
The pipeline is extremely simple to operate. A completely automated delimitation workflow can be initated by typing: 

> python3 HMdelimit.py Imapfile=Imap_file.txt seqfile=seq_file.txt working_dir=Output_Directory

In this case, the pipeline will:

1. Infer a starting phylogeny based on the multiple sequence alignment and Imap file using BPP A01
2. Infer a starting delimitation, and corresponding guide tree using BPP A11
3. Iteratively refine the starting delimitation using the GDI and BPP A00

The output of the pipeline will be:

- A newick tree corresponding to the delimited species
- An Imap file which maps individuals to the delimited species 

### Control of parameters through the command line
Certain parameters of the pipeline can be set directly from the command line. For example, the number of threads used in the BPP calculations can be increased:

> python3 HMdelimit.py Imapfile=Imap_file.txt seqfile=seq_file.txt working_dir=Output_Directory threads=8

### Control of parameters through a Master Control file
The pipeline can also be controlled more directly through a Master Control File (MCF). For example, the mcf could include the following parameters:

> alignment file = D_HL_align.txt
> Imap file = D_HL_imap.txt
> 
> HM guide tree = ((CBC,((NBC,SCA),NCA)),SBC);
> 
> HM mode = merge
> HM GDI threshold = 0.4
> HM decision parameters = any
> 
> thetaprior= 3 0.004 e
> tauprior= 3 0.008
> finetune = 1: 0.01 0.004 0.0003 0.0001 .02
> threads = 2 1 2
> burnin =  5000
> nsample = 8000   
> sampfreq = 4

In this case, the pipeline is provided with a guide tree for the iterative delimitation stage, which enables us to skip the A01 and A11 steps. The bottom 7 parameters of the pipeline are passed directly to the BPP instances. This enables expert users to take granual control of BPP parameters without having to make seperate control files. These parameters will be shared between all modes of BPP, so only mode agnostic parameters such as 'threads' or 'samplefreq' can be used, while mode specific parameters such as 'speciesdelimitation' cannot. As such the pipeline will:

1. Iteratively refine the starting delimitation using the GDI and BPP A00

### Complete control by combining Master Control files, and BPP control files
The individual BPP stages can also be controlled through their own standard BPP control files. In this case, the Master Control file needs to point to the required BPP control files. For example:

> seqfile = align.txt
> Imapfile = imap.txt
> 
> BPP A01 starting phylogeny inference control file = A01.ctl
> BPP A11 starting delimitation control file  = A11.ctl
> BPP A00 HM parameter inference control file = A00.ctl
> 
> threads = 14 2 1
> burnin =   5000
> nsample = 15000   
> sampfreq = 3
> 
> finetune = 1: 5 0.001 0.001 0.001 0.3 0.33 1.0

In such cases, BPP parameters passed from the master control file will be overwritten if the stage specific BPP control file includes the same parameter. However, this still enables us to not have to specify parameters that are shared between instances in each control file, such as 'threads'. 