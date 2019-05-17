#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: python predictTIME.py -I /path2/eVIIS -s strandSpecific -G /path2/refAnnoFile.gtf -B /path2/sample_bam.bam -OF /path2/sample_featureCount.txt -OC /path2/sample_count.txt -OM /path2/sample_fpkm.txt -OT /path2/sample_time.txt
  e.g. python predictTIME.py -I /home/test/eVIIS -s 2 -G /home/test/eVIIS/ref/all.gtf -B /home/test/PT140003/PT140003_mapped_reads.bam -OF /home/test/PT140003/PT140003_featureCount.txt -OC /home/test/PT140003/PT140003_count.txt -OM /home/test/PT140003/PT140003_fpkm.txt -OT /home/test/PT140003/PT140003_time.txt     
       Input:
           -I /path2/eVIIS
           -S 0 (un- stranded), 1 (stranded) and 2 (reversely stranded).   #Indicate if strand-specific read counting should be performed.
           -G /path2/refAnnoFile.gtf
           -B /path2/sample_bam.bam
           
       Output:
           -OF /path2/sample_featureCount.txt
           -OC /path2/sample_count.txt
           -OM /path2/sample_fpkm.txt
           -OT /path2/sample_time.txt
"""


import os 
import sys


def featureCount(sample_bam, strandSpecific, refAnnoFile, sample_featureCount, sample_count):
    cmds = []
    cmds.append("featureCounts -s %d -p -T 6 -a %s -o %s %s" % (strandSpecific, refAnnoFile, sample_featureCount, sample_bam))
    cmds.append("grep -v '^#' %s | cut -f1,6,7 > %s" % (sample_featureCount, sample_count))
    for cmd in cmds:
        os.system(cmd)

def getFPKMtoPredict(refAnnoFile, path2eVIIS, sample_count, sample_fpkm, sample_time):
    cmd = r"python %s/scripts/getFPKMtoPredict.py %s %s %s %s %s" % (path2eVIIS, refAnnoFile, sample_count, sample_fpkm, path2eVIIS, sample_time)
    os.system(cmd)

    

def main():
    path2eVIIS = sys.argv[2]
    strandSpecific = int(sys.argv[4])
    refAnnoFile = sys.argv[6]
    sample_bam = sys.argv[8]
    sample_featureCount = sys.argv[10]
    sample_count = sys.argv[12]
    sample_fpkm = sys.argv[14]
    sample_time = sys.argv[16]
    featureCount(sample_bam, strandSpecific, refAnnoFile, sample_featureCount, sample_count)
    getFPKMtoPredict(refAnnoFile, path2eVIIS, sample_count, sample_fpkm, sample_time)
    
    
if __name__ == "__main__":
    main()