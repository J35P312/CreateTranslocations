import argparse
import subprocess
import sys
import os
sys.path.append(sys.argv[0].replace("wrapper.py",""))
import CreateTranslocations

parser = argparse.ArgumentParser("""This scripts generates a simulated bam file containing structural variants on the input gnome according to the config file""")
parser.add_argument('--fa', type=str, required = True,help="the path to the reference fasta file")
parser.add_argument('--config', type=str, required = True,help="the path to the config file")
parser.add_argument('--prefix', type=str, required = True,help="the prefix of the output files")
parser.add_argument('--insert_size',default="500", type=str,help="library insert size")
parser.add_argument('--read_length',default="150", type=str,help="library read lenght")
parser.add_argument('--insert_std',default="1200", type=str,help="insert size std")
parser.add_argument('--coverage',default=40, type=int,help="mean coverage")
parser.add_argument('--threads',default="1", type=str,help="the number of threads used during the bwa mapping step")
args = parser.parse_args()

#generate the structural variant fa
print(args.insert_size)
simulated_length=CreateTranslocations.main(args)

#run wgsim to generate the fastq files
n=int(round(args.coverage*simulated_length/float(args.read_length)*1/4))
print("simulating reads of haplotype 1")
command=[os.path.join(sys.argv[0].replace("wrapper.py",""),"wgsim/wgsim"),"-d",args.insert_size,"-e","0.0001","-s",args.insert_std,"-N",str(n),"-1",args.read_length,"-2",args.read_length,args.prefix+"_haplotype_1.fa",args.prefix+"_haplotype_1_r1.fq",args.prefix+"_haplotype_1_r2.fq"]
print(" ".join(command))
tmp=subprocess.check_output(command);
print("simulating reads of haplotype 2")
command=[os.path.join(sys.argv[0].replace("wrapper.py",""),"wgsim/wgsim"),"-d",args.insert_size,"-e","0.0001","-s",args.insert_std,"-N",str(n),"-1",args.read_length,"-2",args.read_length,args.prefix+"_haplotype_2.fa",args.prefix+"_haplotype_2_r1.fq",args.prefix+"_haplotype_2_r2.fq"]
print(" ".join(command))
tmp=subprocess.check_output(command);
#append the fastq files
os.system("cat " + args.prefix+"_haplotype_1_r1.fq" + " " + args.prefix+"_haplotype_2_r1.fq" + " > " + args.prefix+"_r1.fq")
os.remove(args.prefix+"_haplotype_1_r1.fq");os.remove(args.prefix+"_haplotype_2_r1.fq")

os.system("cat " + args.prefix+"_haplotype_1_r2.fq" +" " + args.prefix+"_haplotype_2_r2.fq" + " > " + args.prefix+"_r2.fq")
os.remove(args.prefix+"_haplotype_1_r2.fq");os.remove(args.prefix+"_haplotype_2_r2.fq")
#mapping using bwa
command=["bwa mem -t " + args.threads + " -R \'@RG\\tID:CreateTranslocations\\tSM:WGSIM\' "+args.fa+" "+ args.prefix+"_r1.fq " + args.prefix+"_r2.fq > " + args.prefix+".sam"]
print(command)
tmp=subprocess.check_output(command, shell = True);
#use samtools view to convert sam to bam
command=["samtools view -h -S " + args.prefix+".sam | awk  \' /IIIIIIIIIIIIIIIIIIIIIIIIII/ { gsub(\"IIIIII\", \"@@@@@@\"); print $0; next } { print } \' | samtools view -b -h -@ " +args.threads + " -S - > " + args.prefix+".bam" ]
tmp=subprocess.check_output(command,shell = True);
os.remove(args.prefix+".sam")
#samtools sort
command=["samtools","sort","-@",args.threads,args.prefix+".bam",args.prefix+"_sorted"]
tmp=subprocess.check_output(command);
os.remove(args.prefix+".bam")
#bamtools index
command=["samtools","index",args.prefix+"_sorted.bam"]
tmp=subprocess.check_output(command);
