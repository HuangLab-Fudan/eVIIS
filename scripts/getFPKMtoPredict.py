#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: python /path2/eVIIS/scripts/getFPKMtoPredict.py /path2/refAnnoFile /path2/sample_count.txt /path2/sample_fpkm_PC.txt path2eVIIS /path2/sample_time.txt

"""

import sys
import pandas as pd
import numpy as np
import re
from sklearn.externals import joblib


def getPCGeneInformation(refAnnoFile):
    """
    Input:
        refAnnoFile: the reference annotation file used for extracting protein coding genes
    Return:
        pcGeneInfoDict: containing the information only of protein coding genes (key: gene_id, value: gene_name)
    """
    with open(refAnnoFile, 'r') as fobjr:
        allGeneNum = 0
        pcGeneNum = 0
        pcGeneInfoDict = {}
        for eachline in fobjr:
            if eachline[:2] == "##":
                pass  
            else:
                eachline = eachline.strip().split('\t')
                if eachline[2] != "gene":
                    pass
                else:
                    allGeneNum += 1
                    allInfo = eachline[8]
                    allInfo = re.sub('"', '', allInfo)
                    allInfo = allInfo.strip().split(';')
                    itemDict = {}
                    for item in allInfo:
                        if len(item) != 0:
                            item = item.strip().split(' ')
                            itemDict[item[0]] = item[1]
                    if eachline[0] != "chrM" and itemDict['gene_type'] == "protein_coding":
                    ### eliminate chrM genes
                        pcGeneNum += 1
                        pcGeneInfoDict[itemDict['gene_id']] = itemDict['gene_name']               
        print("The total number of genes in %s is: %d " % (refAnnoFile, allGeneNum))
        print("The number of protein_coding genes (excluding chrM genes) in %s is: %d " % (refAnnoFile, pcGeneNum))
        print("There are %d items in pcGeneInfoDict" % (len(pcGeneInfoDict)))
    return pcGeneInfoDict

              
def getFPKM_PC(refAnnoFile, sample_count, sample_fpkm_PC):
    """
    Input:
        sample_count(sample_count.txt): resulted from featureCounts 
    Output:
        sample_fpkm_PC(sample_fpkm_PC.txt): FPKM value (of protein coding gene) transformed from count
    Return:
        fpkmDict: key:Geneid, value:FPKM           
    """
    pcGeneInfoDict = getPCGeneInformation(refAnnoFile)
    df = pd.read_table(sample_count, delim_whitespace = True, header = 0)
    columns = list(df.columns) 
    geneIDs = list(df[columns[0]]) 
    fpkmDict_PC = {}
    fpkmDict_PC['Geneid'] = "FPKM"   
    allCounts_PC = 0
    for pcGeneID in pcGeneInfoDict:
        count = df.iloc[geneIDs.index(pcGeneID), 2]
        allCounts_PC = allCounts_PC + count
    print("the all counts of protein coding genes is: %d" % (allCounts_PC))
    for pcGeneID in pcGeneInfoDict:
        geneIDindex = geneIDs.index(pcGeneID)
        length = df.iloc[geneIDindex, 1]
        count = df.iloc[geneIDindex, 2]
        fpkm = (count * (10**9))/(length * allCounts_PC)
        fpkmDict_PC[pcGeneID] = fpkm          
    fpkmdf_PC = pd.DataFrame(fpkmDict_PC,index=[0])
    fpkmdf_PC = fpkmdf_PC.T
    fpkmdf_PC.to_csv(sample_fpkm_PC, sep= "\t", header = False)
    return fpkmDict_PC


def transStr2Dict(string):
    """
    Input:
        a str format 'dictionary', e.g. "'JAML': '0.77', 'ARHGAP30': '0.27'"
    Return:
        a dict format 'dictionary', e.g. {'JAML': '0.77', 'ARHGAP30': '0.27'}
    """
    string = string.strip().split(',')
    dictionary = {}
    for item in string:
        k = item.strip().split(':')[0].strip()
        k = re.sub("'", '', k)
        v = item.strip().split(':')[1].strip()
        v = re.sub("'", '', v)
        dictionary[k] = v
    return dictionary
 
       
def countLFscore(fpkmDict):
    """
    Input:
        fpkmDict: resulted from getFPKM(sample_count, sample_fpkm_PC)
    Output:
        LFscore   
    """
    lassoSymbolIdDict = "'TNFRSF1B': 'ENSG00000028137', 'PTAFR': 'ENSG00000169403', 'SLAMF8': 'ENSG00000158714', 'ARHGAP30': 'ENSG00000186517', 'CD247': 'ENSG00000198821', 'NCF2': 'ENSG00000116701', 'RASSF5': 'ENSG00000266094', 'CCR4': 'ENSG00000183813', 'TIGIT': 'ENSG00000181847', 'HCLS1': 'ENSG00000180353', 'CD86': 'ENSG00000114013', 'FYB1': 'ENSG00000082074', 'HLA-E': 'ENSG00000204592', 'GPSM3': 'ENSG00000213654', 'HLA-DPB1': 'ENSG00000223865', 'FGD2': 'ENSG00000146192', 'SRGN': 'ENSG00000122862', 'CD5': 'ENSG00000110448', 'JAML': 'ENSG00000160593', 'KLRB1': 'ENSG00000111796', 'SELPLG': 'ENSG00000110876', 'GPR65': 'ENSG00000140030', 'SNX20': 'ENSG00000167208', 'NLRC5': 'ENSG00000140853', 'FMNL1': 'ENSG00000184922', 'LILRB3': 'ENSG00000204577', 'LILRB2': 'ENSG00000131042', 'ZNF831': 'ENSG00000124203', 'GRAP2': 'ENSG00000100351', 'FOXP3': 'ENSG00000049768'"
    lassoSymbolIdDict = transStr2Dict(lassoSymbolIdDict)
    LFDict = "'JAML': '0.77', 'ARHGAP30': '0.27', 'CCR4': '1', 'CD247': '0.62', 'CD5': '0.26', 'CD86': '0.26', 'FGD2': '0', 'FMNL1': '0.25', 'FOXP3': '0.55', 'FYB1': '0.38', 'GPR65': '0.32', 'GPSM3': '0.24', 'GRAP2': '1', 'HCLS1': '0.26', 'HLA-DPB1': '0.25', 'HLA-E': '0.25', 'KLRB1': '0.6', 'LILRB2': '0.41', 'LILRB3': '0.59', 'NCF2': '0.4', 'NLRC5': '0.33', 'PTAFR': '0.38', 'RASSF5': '0.28', 'SELPLG': '0.28', 'SLAMF8': '0.41', 'SNX20': '0.43', 'SRGN': '0.25', 'TIGIT': '0.72', 'TNFRSF1B': '0.25', 'ZNF831': '0.67'"
    LFDict = transStr2Dict(LFDict)
    LFscore = 0
    for symbol in LFDict:
        coe = float(LFDict[symbol])
        fpkm = 0
        geneID = lassoSymbolIdDict[symbol]
        for ensembleID in fpkmDict:
            if geneID not in ensembleID:
                pass
            else:
                fpkm = fpkmDict[ensembleID]
        LFscore += coe * fpkm
    return LFscore

       
def svmPredict(fpkmDict, path2eVIIS):
    """
    Input:
        fpkmDict: resulted from getFPKM(sample_count, sample_fpkm_PC)
        path2eVIIS: the absolute path of eVIIS
    Return:
        svm_pred: the predicted label of this fpkmDict, the value of svm_pred is 0 or 1       
    """
    svmFeaSymbolIdDict = "'CD52': 'ENSG00000169442', 'GBP4': 'ENSG00000162654', 'GBP5': 'ENSG00000154451', 'CD2': 'ENSG00000116824', 'CD8A': 'ENSG00000153563', 'CD8B': 'ENSG00000172116', 'JCHAIN': 'ENSG00000132465', 'CXCL9': 'ENSG00000138755', 'CXCL10': 'ENSG00000169245', 'CXCL13': 'ENSG00000156234', 'GZMK': 'ENSG00000113088', 'GZMA': 'ENSG00000145649', 'CD74': 'ENSG00000019582', 'UBD': 'ENSG00000213886', 'HLA-DRA': 'ENSG00000204287', 'HLA-DRB5': 'ENSG00000198502', 'HLA-DRB1': 'ENSG00000196126', 'HLA-DQA1': 'ENSG00000196735', 'HLA-DQB1': 'ENSG00000179344', 'HLA-DQA2': 'ENSG00000237541', 'HLA-DQB2': 'ENSG00000232629', 'HLA-DOA': 'ENSG00000204252', 'HLA-DPA1': 'ENSG00000231389', 'HLA-DPB1': 'ENSG00000223865', 'IDO1': 'ENSG00000131203', 'PRF1': 'ENSG00000180644', 'CTSW': 'ENSG00000172543', 'CD3E': 'ENSG00000198851', 'CD3D': 'ENSG00000167286', 'CD27': 'ENSG00000139193', 'LAG3': 'ENSG00000089692', 'GZMH': 'ENSG00000100450', 'CCL5': 'ENSG00000271503', 'CEACAM5': 'ENSG00000105388', 'NKG7': 'ENSG00000105374', 'CST7': 'ENSG00000077984', 'IGLL5': 'ENSG00000254709'"    
    svmFeaSymbolIdDict = transStr2Dict(svmFeaSymbolIdDict)
    svmFeaSymsRank = "CXCL9	CCL5	NKG7	CD8A	CXCL10	CXCL13	GZMK	GZMA	HLA-DQA1	CD2	CD3D	CD3E	CST7	CD27	HLA-DPA1	LAG3	GBP5	HLA-DQA2	CD74	HLA-DPB1	HLA-DRB1	CTSW	HLA-DRA	JCHAIN	UBD	HLA-DQB1	HLA-DOA	IDO1	CD52	GBP4	PRF1	HLA-DRB5	HLA-DQB2	GZMH	CD8B	CEACAM5	IGLL5"
    svmFeaSymsRank = svmFeaSymsRank.strip().split('\t') 
    svmInput = []
    for symbol in svmFeaSymsRank:
        fpkm = 0
        geneID = svmFeaSymbolIdDict[symbol]
        for ensembleID in fpkmDict:
            if geneID not in ensembleID:
                pass
            else:
                fpkm = fpkmDict[ensembleID]
        svmInput.append(fpkm)
    svmInput = np.array([svmInput])
    #load scaler
    scaler = joblib.load(path2eVIIS + '/scripts/scaler')
    #scale the svmInput
    svmInput_scale = scaler.transform(svmInput)
    #load model
    svm_model = joblib.load(path2eVIIS + '/scripts/svm_model')
    #use model
    svm_pred = svm_model.predict(svmInput_scale)
    return svm_pred
    

def predict(refAnnoFile, sample_count, sample_fpkm_PC, path2eVIIS, predictFile):
    """
    Input:
        refAnnoFile: the reference annotation file used for extracting protein coding genes
        sample_count(sample_count.txt): resulted from featureCounts
        path2eVIIS: the absolute path of eVIIS
    Output:
        sample_fpkm_PC(sample_fpkm_PC.txt): FPKM value transformed from count
        predictFile: the final result file of eVIIS program
    ### optimal cut-point: LF.Score=81.03
    """
    fpkmDict = getFPKM_PC(refAnnoFile, sample_count, sample_fpkm_PC)
    LFscore = countLFscore(fpkmDict)
    with open(predictFile, 'w') as fobjw:
        print("LF.Score is: %f\n" % (LFscore)) 
        fobjw.write("LF.Score is: %f\n" % (LFscore))
        if LFscore < 81.03:
            print('The LF.Score of this sample is low.\n')
            fobjw.write("The LF.Score of this sample is low.\n")
            print("\nThis sample is predicted as: Immune Exclusion\n")
            fobjw.write("\nThis sample is predicted as: Immune Exclusion\n")
        else:
            print('The LF.Score of this sample is high. This sample will be further classified by SVM model.\n')
            fobjw.write('The LF.Score of this sample is high. This sample will be further classified by SVM model.\n')            
            svm_pred = svmPredict(fpkmDict, path2eVIIS)
            if list(svm_pred)[0] == 0:
                print("\nThis sample is predicted by SVM model as: Immune Stimulation\n")
                fobjw.write("\nThis sample is predicted by SVM model as: Immune Stimulation\n")
            else:
                print("\nThis sample is predicted by SVM model as: Immune Anergy\n")
                fobjw.write("\nThis sample is classified by SVM model as: Immune Anergy\n")

    
def main():
    refAnnoFile = sys.argv[1]
    sample_count = sys.argv[2]
    sample_fpkm_PC = sys.argv[3]
    path2eVIIS = sys.argv[4]
    sample_time = sys.argv[5]
    predict(refAnnoFile, sample_count, sample_fpkm_PC, path2eVIIS, sample_time)
    

if __name__ == "__main__":
    main()
    
    
    

    
        

    

