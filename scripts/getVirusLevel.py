#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: python /path2/eVIIS/scripts/getVirusLevel.py /path2/sample_stringtie.gtf /path2/sample_virus.txt

"""

import re
import sys


def getVirusLevel(assemblyFile, virusLevelFile):
    """
    Input:
        assemblyFile(.gtf): resulted from Stringtie
    Output:
        virusLevelFile(.txt): virus expression level detected in this assemblyFile.gtf 
    """
    virusNames = "chrEBV	CMV	HBV	HCV-1	HCV-2	HIV-1	HIV-2	HPV1	HPV10	HPV100	HPV101	HPV102	HPV103	HPV104	HPV105	HPV106	HPV107	HPV108	HPV109	HPV11	HPV110	HPV111	HPV112	HPV113	HPV114	HPV115	HPV116	HPV117	HPV118	HPV119	HPV12	HPV120	HPV121	HPV122	HPV123	HPV124	HPV125	HPV126	HPV127	HPV128	HPV129	HPV13	HPV130	HPV131	HPV132	HPV133	HPV134	HPV135	HPV136	HPV137	HPV138	HPV139	HPV14	HPV140	HPV141	HPV142	HPV143	HPV144	HPV145	HPV146	HPV147	HPV148	HPV149	HPV15	HPV150	HPV151	HPV152	HPV153	HPV154	HPV155	HPV156	HPV159	HPV16	HPV160	HPV161	HPV162	HPV163	HPV164	HPV165	HPV166	HPV167	HPV168	HPV169	HPV17	HPV170	HPV171	HPV172	HPV173	HPV174	HPV175	HPV178	HPV179	HPV18	HPV180	HPV184	HPV19	HPV197	HPV199	HPV2	HPV20	HPV21	HPV22	HPV23	HPV24	HPV25	HPV26	HPV27	HPV28	HPV29	HPV3	HPV30	HPV31	HPV32	HPV33	HPV34	HPV35	HPV36	HPV37	HPV38	HPV39	HPV4	HPV40	HPV41	HPV42	HPV43	HPV44	HPV45	HPV47	HPV48	HPV49	HPV5	HPV50	HPV51	HPV52	HPV53	HPV54	HPV56	HPV57	HPV58	HPV59	HPV6	HPV60	HPV61	HPV62	HPV63	HPV65	HPV66	HPV67	HPV68	HPV69	HPV7	HPV70	HPV71	HPV72	HPV73	HPV74	HPV75	HPV76	HPV77	HPV78	HPV8	HPV80	HPV81	HPV82	HPV83	HPV84	HPV85	HPV86	HPV87	HPV88	HPV89	HPV9	HPV90	HPV91	HPV92	HPV93	HPV94	HPV95	HPV96	HPV97	HPV98	HPV99	HPV-mCG2	HPV-mCG3	HPV-mCH2	HPV-mFD1	HPV-mFD2	HPV-mFi864	HPV-mFS1	HPV-mKC5	HPV-mKN1	HPV-mKN2	HPV-mKN3	HPV-mL55	HPV-mRTRX7	HPV-mSD2	HTLV-1	KSHV	MCV	SV40	HPV157	HPV158	HPV187	HPV200	HCV-3	HCV-4	HCV-5	HCV-6	Hepatitis delta virus	HEV	EBV-1	EBV-2	Human herpesvirus 1 (CMV-1)	Human herpesvirus 2  (CMV-2)	Human herpesvirus 3  (CMV-3)	SV41	HTLV-2	HTLV-3"
    virusNames = virusNames.split('\t')
    virusTPMdict = {}
    for virus in virusNames:
        virusTPMdict[virus] = []
    with open(assemblyFile, 'r') as fobjr, open(virusLevelFile, 'w') as fobjw:
        for eachline in fobjr:
            eachlineList = eachline.strip().split('\t')
            chromNum = eachlineList[0]
            if chromNum not in virusNames:
                pass
            else:     
                infoList = eachlineList[8].strip().split(';')          
                if len(infoList) != 6:
                    pass
                else:
                    TPM = infoList[4]
                    TPM = re.sub('"', '', TPM)
                    virusTPMdict[chromNum].append(float(TPM.strip().split(' ')[1]))
                    fobjw.write(eachline)
        fobjw.write('\n')
        fobjw.write('\n')
        fobjw.write("-------------------- Virus expression level in this sample is: --------------------\n")
        print("-------------------- Virus expression level in this sample is: --------------------\n")
        fobjw.write("Virus\tTPM\t\n")
        num = 0
        for virus in list(virusTPMdict.keys()):
            TPMs = virusTPMdict[virus]
            if len(TPMs) != 0:
                fobjw.write("%s\t%f\n" % (virus, max(TPMs)))
                print("%s\t%f\n" % (virus, max(TPMs)))
            else:
                num += 1
        if num == len(virusTPMdict):
            fobjw.write("No virus was detected in this sample.")
            print("No virus was detected in this sample.")
    
        else:
            pass
            
def main():
    sample_stringtie = sys.argv[1]
    sample_virus = sys.argv[2]
    getVirusLevel(sample_stringtie, sample_virus)
    
    
                
if __name__ == "__main__":
    main()   

        
    
    
                
        
        