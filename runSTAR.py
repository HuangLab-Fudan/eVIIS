#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

### Create STAR index
Usage: python runSTAR.py -I /path2/eVIIS -index -SI_dir /path2/star_index -G /path2/refAnnoFile.gtf -GA /path2/refFastaFile.fasta
  e.g. python runSTAR.py -I /home/test/eVIIS -index -SI_dir /home/test/eVIIS/star_index -G /home/test/eVIIS/ref/all.gtf -GA /home/test/eVIIS/ref/all.fasta   
       Input: 
           -I /path2/eVIIS
           -index
           -G /path2/refAnnoFile.gtf
           -GA /path2/refFastaFile.fasta   
       Output:
           -SI_dir /path2/star_index (files in SI_dir)
           
### STAR to map RNA-seq data
Usage: python runSTAR.py -I /path2/eVIIS -pass -SI_dir /path2/star_index -f1 seq1.fq.gz -f2 seq2.fq.gz -fq_dir /path2/seqs -GA /path2/refFastaFile.fasta -O /path2/output -S sample
  e.g. python runSTAR.py -I /home/test/eVIIS -pass -SI_dir /home/test/eVIIS/star_index -f1 PT140003_R1.fq.gz -f2 PT140003_R2.fq.gz -fq_dir /home/D/HNSC/fa -GA /home/test/eVIIS/ref/all.fasta -O /home/test -S PT140003
       Input:
           -I /path2/eVIIS
           -pass
           -SI_dir /path2/star_index
           -f1 seq1.fq.gz
           -f2 seq2.fq.gz
           -fq_dir /path2/seqs
           -GA /path2/refFastaFile.fasta
           -S sample  
       Output:
           -O /path2/output/ (files in this path, including sample_mapped_reads.bam)
"""

import os 
import sys


def makeDir(path):
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
  

def starIndex(path2eVIIS, path2index, refAnnoFile, refFastaFile):
    cmd= r"sh %s/scripts/star_index.sh %s %s %s" % (path2eVIIS, path2index, refFastaFile, refAnnoFile) 
    os.system(cmd)
    

def starAlign(path2eVIIS, path2index, fqgz1, fqgz2, path2fqgz, refFastaFile, path2output, sampleName):
    cmd = r"sh %s/scripts/two_pass.sh %s %s %s %s %s %s %s" % (path2eVIIS, path2index, refFastaFile, sampleName, path2fqgz, path2output, fqgz1, fqgz2)
    os.system(cmd)
    
    
def main():
    if len(sys.argv) == 10:
        path2eVIIS = sys.argv[2]
        path2index = sys.argv[5]
        refAnnoFile = sys.argv[7]
        refFastaFile = sys.argv[9]
        makeDir(path2index)
        print("The STAR genome index will be created using this command.\n")
        starIndex(path2eVIIS, path2index, refAnnoFile, refFastaFile)
    elif len(sys.argv) == 18:
        path2eVIIS = sys.argv[2]
        path2index = sys.argv[5]
        fqgz1 = sys.argv[7]
        fqgz2 = sys.argv[9]
        path2fqgz = sys.argv[11]
        refFastaFile = sys.argv[13]
        path2output = sys.argv[15]
        sampleName = sys.argv[17]
        makeDir(path2output)
        print("Tow pass of STAR will be executed using this command.\n")
        starAlign(path2eVIIS, path2index, fqgz1, fqgz2, path2fqgz, refFastaFile, path2output, sampleName)
    else:
        print("Some variables may not be defined, please check it carefully.\n")
        
    
    
if __name__ == "__main__":
    main()