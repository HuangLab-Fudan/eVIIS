#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: python detectVirus.py -I /path2/eVIIS -G /path2/refAnnoFile.gtf -B /path2/sample_bam.bam -OS /path2/sample_stringtie.gtf -OV /path2/sample_virus.txt
  e.g. python detectVirus.py -I /home/test/eVIIS -G /home/test/eVIIS/ref/all.gtf -B /home/test/PT140003/PT140003_mapped_reads.bam -OS /home/test/PT140003/PT140003_stringtie.gtf -OV /home/test/PT140003/PT140003_virus.txt     
       Input:
           -I /path2/eVIIS
           -G /path2/refAnnoFile.gtf
           -B /path2/sample_bam.bam
           
       Output:
           -OS /path2/sample_stringtie.gtf
           -OV /path2/sample_virus.txt    
"""

import os 
import sys


def stringtieAssemble(sample_bam, refAnnoFile, sample_stringtie):
    ###-f, -p
    cmd = r"stringtie %s -f 0.1 -o %s -p 15 -G %s" % (sample_bam, sample_stringtie, refAnnoFile)
    os.system(cmd)
    
    
def getVirusLevel(path2eVIIS, sample_stringtie, sample_virus):
    cmd = r"python %s/scripts/getVirusLevel.py %s %s" % (path2eVIIS, sample_stringtie, sample_virus)
    os.system(cmd)
    
def main():
    path2eVIIS = sys.argv[2]
    refAnnoFile = sys.argv[4]
    sample_bam = sys.argv[6]
    sample_stringtie = sys.argv[8]
    sample_virus = sys.argv[10]
    stringtieAssemble(sample_bam, refAnnoFile, sample_stringtie)
    getVirusLevel(path2eVIIS, sample_stringtie, sample_virus)
    
    
if __name__ == "__main__":
    main()
    
    
        
