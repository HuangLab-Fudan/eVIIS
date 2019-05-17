#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
### Create STAR index
Usage: python eVIIS.py -I /path2/eVIIS -index -SI_dir /path2/star_index -G /path2/refAnnoFile.gtf -GA /path2/refFastaFile.fasta
       Input: 
           -I /path2/eVIIS
           -index
           -G /path2/refAnnoFile.gtf
           -GA /path2/refFastaFile.fasta   
       Output:
           -SI_dir /path2/star_index (files in SI_dir)
           
### Map and assemble RNA-seq data, and estimate viral expression level, predict immune status (TIME)
Usage: python eVIIS.py -I /path2/eVIIS -pass -SI_dir /path2/star_index -f1 seq1.fq.gz -f2 seq2.fq.gz -fq_dir /path2/seqs -G /path2/refAnnoFile.gtf -GA /path2/refFastaFile.fasta -O /path2/output -S sample -s strandSpecific
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
           -O /path2/output (files in this path, including:/sample/
                              sample_mapped_reads.bam, 
                              sample_stringtie.gtf
                              sample_virus.txt
                              sample_featureCount.txt
                              sample_count.txt
                              sample_fpkm.txt
                              sample_time.txt)
"""

import os
import sys
import runSTAR as rs
import detectVirus as dv
import predictTIME as pt


def makeDir(path):
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
        

def main():
    if len(sys.argv) == 10:
        path2eVIIS = sys.argv[2]
        path2index = sys.argv[5]
        refAnnoFile = sys.argv[7]
        refFastaFile = sys.argv[9]
        makeDir(path2index)
        print("The STAR genome index will be created using this command.\n")
        rs.starIndex(path2eVIIS, path2index, refAnnoFile, refFastaFile)
    elif len(sys.argv) == 22:
        path2eVIIS = sys.argv[2]
        path2index = sys.argv[5]
        fqgz1 = sys.argv[7]
        fqgz2 = sys.argv[9]
        path2fqgz = sys.argv[11]
        refAnnoFile = sys.argv[13]
        refFastaFile = sys.argv[15]
        path2output = sys.argv[17]
        sampleName = sys.argv[19]
        strandSpecific = int(sys.argv[21])
        makeDir(path2output)
        print("Tow pass of STAR will be executed using this command.\n")
        rs.starAlign(path2eVIIS, path2index, fqgz1, fqgz2, path2fqgz, refFastaFile, path2output, sampleName)
    
        sample_bam = path2output+'/'+sampleName+'/'+sampleName+'_mapped_reads.bam'
        sample_stringtie = path2output+'/'+sampleName+'/'+sampleName+'_stringtie.gtf'
        sample_virus = path2output+'/'+sampleName+'/'+sampleName+'_virus.txt'
        sample_featureCount = path2output+'/'+sampleName+'/'+sampleName+'_featureCount.txt'
        sample_count = path2output+'/'+sampleName+'/'+sampleName+'_count.txt'
        sample_fpkm = path2output+'/'+sampleName+'/'+sampleName+'_fpkm.txt'
        sample_time = path2output+'/'+sampleName+'/'+sampleName+'_time.txt'
    
        dv.stringtieAssemble(sample_bam, refAnnoFile, sample_stringtie)
        dv.getVirusLevel(path2eVIIS, sample_stringtie, sample_virus)
        
        pt.featureCount(sample_bam, strandSpecific, refAnnoFile, sample_featureCount, sample_count)
        pt.getFPKMtoPredict(refAnnoFile, path2eVIIS, sample_count, sample_fpkm, sample_time)

    else:
        print("Some variables may not be defined, please check it carefully.\n")

if __name__ == "__main__":
    main()

