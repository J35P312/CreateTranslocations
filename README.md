# CreateTranslocations
1. use the install script to download wgsim and SURVIVOR, make sure to have bwa,samtools, and sambamba and installed
2. run SURVIVOR to create a config file : ./SURVIVOR/Debug/SURVIVOR simSV config.txt
3. edit the config.txt file to suit your needs
4. run wrapper.py

# dependencies

    bwa, sambamba and samtools

# Command line

    To produce a help message, type

        python wrapper.py --help


# note

Crate translocation will create a diploid genome, hence the analysis is run twice, generating half the set coverage each run.
