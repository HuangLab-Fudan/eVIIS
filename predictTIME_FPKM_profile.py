#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

### Create STAR index
Usage: python predictTIME_FPKM_profile.py -I /path2/eVIIS -fp /path2/FPKM_profile.txt -OT_dir /path2/sample_time -OR /path2/result_all.txt
       Input: 
           -I /path2/eVIIS
           -fp /path2/FPKM_profile.txt
       Output:
           -OT_dir /path2/sample_time
           -OR /path2/result_all.txt
Note:
     FPKM_profile.txt:  Each row of FPKM_profile.txt is a gene, each column of FPKM_profile.txt is a sample. The first row is sample names, the first column is Geneids.

"""         



import sys
import pandas as pd
import numpy as np
import re
from sklearn.externals import joblib


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
        fpkmDict: 
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
        fpkmDict: 
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
    

def predictSingleSample(fpkmDict, path2eVIIS, predictFile):
    """
    Input:
        fpkmDict
        path2eVIIS: the absolute path of eVIIS
    Output:
        predictFile: the final result file of TIME prediction
    ### optimal cut-point: LF.Score=81.03
    """
    fpkmDict = fpkmDict
    LFscore = countLFscore(fpkmDict)
    with open(predictFile, 'w') as fobjw:
        print("LFscore is: %f\n" % (LFscore)) 
        fobjw.write("LFscore is: %f\n" % (LFscore))
        if LFscore < 81.03:
            print('The LFscore of this sample is low.\n')
            fobjw.write("The LFscore of this sample is low.\n")
            print("\nThis sample is predicted as: Immune Exclusion\n")
            fobjw.write("\nThis sample is predicted as: Immune Exclusion\n")
        else:
            print('The LFscore of this sample is high. This sample will be further classified by SVM model.\n')
            fobjw.write('The LFscore of this sample is high. This sample will be further classified by SVM model.\n')            
            svm_pred = svmPredict(fpkmDict, path2eVIIS)
            if list(svm_pred)[0] == 0:
                print("\nThis sample is predicted by SVM model as: Immune Stimulation\n")
                fobjw.write("\nThis sample is predicted by SVM model as: Immune Stimulation\n")
            else:
                print("\nThis sample is predicted by SVM model as: Immune Anergy\n")
                fobjw.write("\nThis sample is classified by SVM model as: Immune Anergy\n")


def predictFPKMprofile(path2eVIIS, fpkm_profile, sample_time_filepath, resultFile):  
    """
    Input:
        path2eVIIS: the absolute path of eVIIS, to get the svm_model
        fpkm_profile (FPKM_profile.txt):  Each row of fpkm_profile is a gene, each column of fpkm_profile is a sample. The first row is sample names, the first column is Geneids.
    Output:
        sample_time_filepath: the path of TIME prediction result files
        resultFile (result_all.txt): the integration file of the prediction results of all samples in the fpkm_profile.
    """
    fpkm_df = pd.read_table(fpkm_profile, delim_whitespace = True, header=0, index_col=0)
    print(fpkm_df.shape)
    sample_names = list(fpkm_df.columns)
    gene_ids = list(fpkm_df.index)
    for sample_name in sample_names:
        sample_time = sample_time_filepath+'/'+sample_name+"_time.txt"
        fpkmDict = {}
        for gene_id in gene_ids:
            fpkmDict[gene_id] = fpkm_df.loc[gene_id, sample_name]
        print(len(fpkmDict))
        predictSingleSample(fpkmDict, path2eVIIS, sample_time)
    ###integrate the result of each sample into a file "result_all.txt"
    with open(resultFile, 'w') as fobjw:
        fobjw.write("Sample\tLF.Score\tTIME\n")
        for sample_name in sample_names:
            sample_time = sample_time_filepath+'/'+sample_name+"_time.txt"
            predict_results = []
            with open(sample_time, 'r') as fobjr:
                for eachline in fobjr:
                    eachline= eachline.strip().split(":")
                    if len(eachline) == 2:
                        predict_results.append(eachline[1])
                    else:
                        pass
            newline = sample_name + '\t' + '\t'.join(predict_results) + '\n'
            fobjw.write(newline) 
            
            
def main():
    #-I
    path2eVIIS = sys.argv[2]
    #-fp
    fpkm_profile = sys.argv[4]
    #-OT_dir
    sample_time_filepath = sys.argv[6]
    #-OR
    resultFile = sys.argv[8]
    predictFPKMprofile(path2eVIIS, fpkm_profile, sample_time_filepath, resultFile)


if __name__ == "__main__":
    main()
        
   






