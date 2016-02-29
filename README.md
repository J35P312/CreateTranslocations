# CreateTranslocations
1. use the install script to download wgsim, make sure to have bwa,samtools and bamtools installed
2. edit the config.txt file to suit your needs
3. python wrapper.py

the config file, an exmple

the chromosome line is used to set which chromosomes that are to be simulated, all variants are spred randomly accross all selected chromosomes

chromosomes=chr18,chr20

the fraction of homozygous events(diploi genome)
homozygozity_rate=0.5

the minimum distance between events on the same haplotype
min_dist=100000

translocations, these are insertions of sequence from a randomly chosen chromosome, the first number is the number of each kind of translocation to simulate, the second is the smallest inserted sequence, and the third the largest. The 4th and last number is the inversion frequnce, i.e ho often the inverted sequence in inverted

balanced_inter_chromosomal_translocations=300,2000,100000,0.5

unbalanced_inter_chromosomal_translocations=300,2000,100000,0.5

number of variants, minimum size, and maximum size

inversions=100,2000,100000

deletions=300,2000,100000

tandem_duplications=100,2000,100000

modify the numbers to suit your needs

within the wrapper.py script you can set parameters such as coverage and insert size.
The output of the pipeline is a sorted and indexed bam file, as well as a tab file describing each variant
