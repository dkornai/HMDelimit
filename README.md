# HMDelimit
This repo holds the python source files for a species delimitation pipeline based on the work of Leache et al 2019.
# Methods:
The pipeline performs multilocus species delimitation by combining the use of BPP, and a metric of speciation 
called the GDI. 
# Examples:
The pipeline can be extremely simple to operate. 
A completely automated delimitation workflow can be initated by typing: 

python3 HMdelimit.py D_ROT_align.txt D_ROT_imap.txt Output_directory

The pipeline can also be operated more granularly with the help of a Master control file.