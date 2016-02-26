import sys
import os
import re
import random
import operator
import copy
import argparse

def invert_sequence(sequence):
    inverted_sequence=sequence[::-1]
    return(inverted_sequence)

def check_overlap(structural_variant,haplotype,min_distance):
    overlap=False
    for event in structural_variant:
        for variant in haplotype[event[0]]:
            if variant[1] == event[1]:
                #the added variant must be positioned min_distance or greater from all variants within the haplotype
                if (abs(variant[2]-event[3]) < min_distance) and (abs(variant[3]-event[2]) < min_distance ):
                    overlap=True

                #the added variant should not overlap any variant within the haplotype
                if variant[3] >= event[2] and event[3] >= variant[2]:
                    overlap = True
    return(overlap)

#generate a list containing the variants
def generate_variant_lists(data,chromosome_length):
    n_chomosomes=len(data["chromosomes"])
    number_of_variants =0;
    number_of_variants += data["balanced_inter_chromosomal_translocations"][0]
    number_of_variants += data["unbalanced_inter_chromosomal_translocations"][0]
    number_of_variants += data["inversions"][0]
    number_of_variants += data["deletions"][0]
    number_of_variants += data["tandem_duplications"][0]
    
    variants=["balanced_inter_chromosomal_translocations","unbalanced_inter_chromosomal_translocations","inversions","deletions","tandem_duplications"]
    genome=[{},{}]
    collaped_haplotypes={}
    for chromosome in data["chromosomes"]:
        genome[0][chromosome]=[]
        genome[1][chromosome]=[]
        collaped_haplotypes[chromosome]=[]
    while number_of_variants > 0:
        #decide if the variant is homozygous
        homozygous = 0
        if(random.random() <= float(data["homozygozity_rate"])):
            homozygous = 1
        else:
            #if the variant is not homoyzygous, then choose which haplotype to affect
            haplotype=random.randint(0, 1)
            
        #choose a random chromsome
        chromosome=random.choice(data["chromosomes"])
        #select variant type
        variant=""
        while variant == "":
            variant_type=random.choice(variants)
            if(data[variant_type][0] > 0):
                variant=variant_type;
        #select the start and end position of the variant
        if variant_type is "balanced_inter_chromosomal_translocations" or variant_type is "unbalanced_inter_chromosomal_translocations":
            invert=0
            if(random.random() < data[variant_type][-1]):
                invert=1
            #first select the position of chromosomeA
            startA=random.randint(0, chromosome_length[chromosome])           
            #then select the region on chromosome B
            chromosomeB=chromosome
            while chromosomeB == chromosome:
                chromosomeB=random.choice(data["chromosomes"])
            startB=random.randint(0, chromosome_length[chromosomeB])
            length=random.randint(data[variant_type][1], data[variant_type][2])
            endB=startB+length-1
            structural_variant=[[chromosome,chromosomeB,startA,startB,[startA,endB,invert],"BND"]]
            #if the variant is a balanced translocation, a deletion is generated on chromosome B       
            if variant_type == "balanced_inter_chromosomal_translocations":
                structural_variant += [[chromosomeB,chromosomeB,startB,endB,"1","DEL"]]
 
            #the FT database stores variants according to lexiographic order
            if(chromosome > chromosomeB):
                chromosome=chromosomeB
        else:
            posA=random.randint(0, chromosome_length[chromosome])
            length=random.randint(data[variant_type][1], data[variant_type][2])
            var="DEL"
            if variant_type is "inversions":
                var ="INV"
            elif variant_type is "tandem_duplications":
                var="TDUP"
            structural_variant=[[chromosome,chromosome,posA,posA+length,"1",var]]
        #make sure that the variant do not overlap any other variant
        if homozygous:
            overlap=check_overlap(structural_variant,genome[0],data["min_distance"])
            if not overlap:
                overlap=check_overlap(structural_variant,genome[1],data["min_distance"])
        else:
            overlap=check_overlap(structural_variant,genome[haplotype],data["min_distance"])
            
            
        #if the variant do not overlap, then add it to the genome
        if not overlap:
            number_of_variants += -1
            data[variant_type][0] += -1
            if(homozygous):
                for event in structural_variant:
                    genome[0][event[0]]+=[event]
                    genome[1][event[0]]+=[event]
            else:
                for event in structural_variant:
                    genome[haplotype][event[0]]+=[event]
            for event in structural_variant:
                collaped_haplotypes[event[0]] += [event]

    return(genome,collaped_haplotypes)

def main(args):
    #argument 1 is the fasta file, argument 2 is the variant configfile argument3 is the prefix of the output   bam file
    fafile=args.fa
    bedfile=args.config
    prefix=args.prefix
    variants=[]
    #add all the variants to the variants list
    data={}
    homozygozity_rate=0
    for line in open(bedfile):
        content=line.strip().split("=")
        if content[0] == "chromosomes":
            data["chromosomes"]=content[1].split(",")
        elif content[0] == "min_dist":
            data["min_distance"]=int(content[1])
        elif content[0] == "homozygozity_rate":
            data["homozygozity_rate"]=float(content[1])
        elif content[0] == "balanced_inter_chromosomal_translocations":
            data["balanced_inter_chromosomal_translocations"]= [ float(x) for x in content[1].split(",") ]
        elif content[0] == "unbalanced_inter_chromosomal_translocations":
            data["unbalanced_inter_chromosomal_translocations"]= [ float(x) for x in content[1].split(",") ]
        elif content[0] == "inversions":
            data["inversions"]= [ int(x) for x in content[1].split(",") ]
        elif content[0] == "deletions":
            data["deletions"]= [ int(x) for x in content[1].split(",") ]
        elif content[0] == "tandem_duplications":
            data["tandem_duplications"]= [ int(x) for x in content[1].split(",") ]

    sequence={}
    #read the fast file
    with open(fafile, 'r+') as f:
        reference = f.read()
    split_reference=reference.split(">")
    del reference
    del split_reference[0]
    #store the reference as a dictionary
    chromosome_len={}
    chromosomes=[]
    simulated_bases=0
    for chromosome in split_reference:
        content=chromosome.split("\n",1)
        if content[0] in data["chromosomes"]:
            #print(content)
            sequence[content[0]]=content[1].replace("\n","")
            chromosome_len[content[0]]=len(sequence[content[0]])
            simulated_bases += len(sequence[content[0]])
    del split_reference
    #generate the variants
    haplotypes,collapsed_variants=generate_variant_lists(data,chromosome_len)
    #print the generated variants to prefix.db
    with open(prefix+'.db', 'w') as output_db:
        for chromosome in sorted(collapsed_variants):
            for variants in sorted(collapsed_variants[chromosome], key=operator.itemgetter(2)):
                #the db events must be written in lexiographic order(otherwise BND event will cause problems)
                if variants[0] <= variants[1]:
                    #each translocation have 2 breakpoints
                    if(variants[5] == "BND"):
                        output_db.write( "\t".join([variants[0],variants[1],str(variants[2]),str(variants[2]),str(variants[3]),str(variants[3]),variants[5],"1"])+"\n")
                        output_db.write( "\t".join([variants[0],variants[1],str(variants[2]),str(variants[2]),str(variants[4][1]),str(variants[4][1]),variants[5],"1"])+"\n")
                    else:
                        output_db.write( "\t".join([variants[0],variants[1],str(variants[2]),str(variants[3]),str(variants[2]),str(variants[3]),variants[5],"1"])+"\n")
                else:
                    if(variants[5] == "BND"):
                        output_db.write( "\t".join([variants[1],variants[0],str(variants[3]),str(variants[3]),str(variants[2]),str(variants[2]),variants[5],"1"])+"\n")
                        output_db.write( "\t".join([variants[1],variants[0],str(variants[4][1]),str(variants[4][1]),str(variants[2]),str(variants[2]),"1",variants[5],"1"])+"\n")
                    
                a=1
    del collapsed_variants

    #create one fasta file per haplotype
    generated_genome=[]
    i=1
    for haplotype in haplotypes:
        with open(prefix+"_haplotype_{0}.fa".format(i), 'w') as output_fa:
            i +=1
            genome=""
            for chromosome in sorted(haplotype):
                print chromosome
                previous_bnd=0
                output_fa.write(">"+chromosome+"\n")
                for variant in sorted(haplotype[chromosome], key=operator.itemgetter(2)):
                    start=variant[2];end=variant[3]
                    genome += sequence[chromosome][previous_bnd:start]
                    if variant[-1] == "DEL":
                        #deletion, skip the sequence between start and stop
                        pass
                    elif variant[-1] == "INV":
                        #inversion, add the inverted sequnce between start and stop to the genom
                        genome += invert_sequence(sequence[chromosome][start:end])
                    elif variant[-1] == "TDUP":
                        #tdup, add an extra copy of start-stop directly after the first copy
                        genome += sequence[chromosome][start:end]
                        genome += sequence[chromosome][start:end]
                    elif variant[-1] == "BND":
                        #add a segment of chromosome b to chromosome a, the segment may be inverted
                        if not variant[4][-1]:              
                            genome += sequence[variant[1]][variant[3]:variant[4][1]]
                        else:
                            genome +=invert_sequence(sequence[variant[1]][variant[3]:variant[4][1]])
                        end=start
                    if(len(genome) > 500):                
                        output_fa.write(genome+"\n")
                        genome=""
                    previous_bnd=end
            if previous_bnd < chromosome_len[chromosome]:
                genome += sequence[chromosome][previous_bnd:]
            output_fa.write(genome+"\n")
    return(simulated_bases)    
#parser = argparse.ArgumentParser("""This scripts generates structural variants on the input genome according to the config file""")
#parser.add_argument('--fa', type=str, required = True,help="the path to the reference fasta file")
#parser.add_argument('--config', type=str, required = True,help="the path to the config file")
#parser.add_argument('--prefix', type=str, required = True,help="the prefix of the output files")
#args = parser.parse_args()
#main(args)
