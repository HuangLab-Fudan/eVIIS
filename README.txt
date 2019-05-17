
eVIIS: a program for evaluating Viral Infection and Immune Status
-----------------------------------------------------------------
eVIIS manuscript now published in XXX. Open Access: www.XXX


What is the eVIIS
-----------------
Human oncogenic viruses, identified by a variety of traditional virological techniques, have been nominated causing 10–15% of human cancers worldwide. It’s important to estimate the virus infection status of tumor patients. Only a limited number of programs are available for detecting virus expression level or even for assessing immune status based on RNA-seq data.
Here, we construct the eVIIS, a software package that can estimate the viral infection and immune status from high-throughput RNA sequencing (RNA-seq) data. eVIIS processes assembled transcripts from the STAR aligner and StringTie assembler to detect the virus expression level in corresponding sample, processes quantification information of genes from featureCounts to predict the immune status of corresponding sample. eVIIS provides the stepwise or one-step analysis of RNA-seq data for one sample. Users can choose the appropriate mode of eVIIS according to different data formats or analysis purposes.  


Implementation and Dependencies
-------------------------------
eVIIS was developed with python (3.6.0) and shell (bash) language. Before running the program, it is necessary to check or install python packages as follows:
*os
*sys
*res
*pandas
*numpy
*sklearn

Moreover, eVIIS works based on the STAR, StringTie and featureCounts, so these tools also should be installed and their pathway should be added in ~/.bashrc
* STAR (version <= 2.5)
* StringTie (version <= 1.2.3)
* featureCounts (version >= 1.5.0)


eVIIS Installation
------------------
Download the package and then unzip it in Linux (CentOS or Ubuntu) 


Documentation
-------------
This chapter provides detailed commands arguments and description of output. These commands are labeled with ‘Usage’.

Files Needed:
------------
1.	Human (hg38.fasta) and virus genome sequences and GTF Files (we recommend GENCODE, and the program will report an error if you use GTF from UCSC) are used to generate STAR index.
        (We have provided two reference files named as "all.fasta" and "all.gtf" in the "ref" folder, which are consisted of genome sequences and annotation information for human and more than 200 types of virus. The names and GenBank IDs of these viruses are listed in the file “virus_list.txt”. Users can apply these two files to generate STAR index directly.)
2.	Raw data of RNA-seq (fasta.gz)


Commands and arguments:
----------------------
*	Note: the absolute path is necessary to perform scripts 
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
################################################################### Analyze RNA-seq data of one sample using eVIIS step-by-step ###################################################################
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Step 1, Mapping the RNA-seq data.
Usage: python runSTAR.py [OPTIONS] 

The arguments of runSTAR.py are as followings and if you want to work with single-read, please see more detail of STAR on https://github.com/alexdobin/STAR:
	-f1 <FASTA1>  
               Using Illumina paired-end reads, and the name of read1 has to be supplied.
	-f2 <FASTA2>  
               Using Illumina paired-end reads, and the name of read2 has to be supplied.
	-fq_dir <fastq dir> 
               Specifies path to files containing the sequences to be mapped.
	-G <path_to_gtf file> 
                Specifies the path to the file with annotated transcripts in the standard GTF format.
	-GA <path_to_fasta files> 
                Specified one or more FASTA files with the genome reference sequences.
	-O <output dir>   
                Specifies path to the directory of result files.
	-pass 
                Running STAR in the 2-pass alignment mode.
	-index
                Generating genome index of STAR with default settings.
	-SI_dir < genome index dir>  
                specifies path to the genome directory where genome indexes where generated.
	-I <path>  
                Specifies path to the directory where the eVIIS installation.
	-S <sample>   
                Name of sample

### Generating genome indexes ###
        Usage: python runSTAR.py -I /path2/eVIIS -index -SI_dir /path2/star_index -G /path2/refAnnoFile.gtf -GA /path2/refFastaFile.fasta
### Running STAR in the 2-pass mode [Kahles et al., 2018, Cancer Cell 34, 1–14] ###
        Usage: python runSTAR.py -I /path2/eVIIS -pass -SI_dir /path2/star_index -f1 seq1.fq.gz -f2 seq2.fq.gz -fq_dir /path2/seqs -GA /path2/refFastaFile.fasta -O /path2/output -S sample

“””
Files Needed:
             refFastaFile.fasta
             refAnnoFile.gtf 
             seq1.fq.gz
             seq2.fq.gz
Output: 
       sample_mapped_reads.bam   (***)
“””

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Step 2, Detecting viral infection.
Usage: python detectVirus.py [OPTIONS]
In this step, StringTie is firstly employed to assemble the alignment BAM file. (Users can access https://ccb.jhu.edu/software/stringtie to see more details about StringTie.) We then quantify viral expression level in one sample based on the TPM values in the assembly file. 

The arguments of detectVirus.py are as followings:
        -B <path_to_bam file>
                Specifies the path to the bam file resulted from STAR aligning. This file is the input of StringTie.
	-G <path_to_gtf file> 
                Specifies the path to the file with annotated transcripts in the standard GTF format.
	-I <path>  
                Specifies path to the directory where the eVIIS installation.
	-OS <path_to_stringtie file>   
                Specifies the path and name for the result file of StringTie assembly analysis. 
        -OV <path_to_virus file>
                (***) Specifies the path and name for the result file of viral infection detection. The detected viruses and expression levels (TPM values) are listed in this file.

### Evaluating virus expression level ###
        Usage: python detectVirus.py -I /path2/eVIIS -G /path2/refAnnoFile.gtf -B /path2/sample_mapped_reads.bam -OS /path2/sample_stringtie.gtf -OV /path2/sample_virus.txt

“””
Files Needed:
             refAnnoFile.gtf
             sample_mapped_reads.bam
Output: 
       sample_stringtie.gtf
       sample_virus.txt   (***)
“””

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Step 3, Predicting immune status.
Usage: python predictTIME.py [OPTIONS]
In this step, the read counts of gene level can be calculated by featureCounts (Please see http://subread.sourceforge.net/). Then, the count number of each gene will be transformed into FPKM value based on a formulas (Please refer to The GDC mRNA quantification analysis pipeline, https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm). 

The arguments of predictTIME.py are as followings:
        -B <path_to_bam file>
                Specifies the path to the bam file resulted from STAR aligning. This file is the input of featureCounts.
        -s <int>
                A parameter for featureCounts to indicate if strand-specific read counting should be performed. Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).  
	-G <path_to_gtf file> 
                Specifies the path to the file with annotated transcripts in the standard GTF format.
	-I <path>  
                Specifies path to the directory where the eVIIS installation.
	-OF <path_to_featureCount file>   
                Specifies the path and name for the result file of featureCounts, the read counts and annotation information are included in this file.
        -OC <path_to_count file>
                Specifies the path and name for count file which is derived from featureCount file and only contains geneID and count number of each gene.  
	-OM <path_to_fpkm file>   
                Specifies the path and name for the FPKM file.
        -OT <path_to_time file>
                (***) Specifies the path and name for the final result file of predicting immune status. The calculated LF.Score and predicted TIME are contained in this file.
                
### Counting LS.Score and predicting TIME  ###
        Usage: python predictTIME.py -I /path2/eVIIS -s 0 -G /path2/refAnnoFile.gtf -B /path2/sample_mapped_reads.bam -OF /path2/sample_featureCount.txt -OC /path2/sample_count.txt -OM /path2/sample_fpkm.txt -OT /path2/sample_time.txt

“””
Files Needed:
             refAnnoFile.gtf
             sample_mapped_reads.bam
Output: 
       sample_featureCount.txt
       sample_count.txt
       sample_fpkm.txt
       sample_time.txt   (***)
“””

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#################################################################### Analyze RNA-seq data of one sample using eVIIS by one step ###################################################################
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
One step, Mapping the RNA-seq data, Detecting viral infection, and Predicting immune status.
Usage: python eVIIS.py [OPTIONS] 

The arguments of eVIIS.py are as followings:
	-f1 <FASTA1>  
               Using Illumina paired-end reads, and the name of read1 has to be supplied.
	-f2 <FASTA2>  
               Using Illumina paired-end reads, and the name of read2 has to be supplied.
	-fq_dir <fastq dir> 
               Specifies path to files containing the sequences to be mapped.
	-G <path_to_gtf file> 
                Specifies the path to the file with annotated transcripts in the standard GTF format.
	-GA <path_to_fasta files> 
                Specified one or more FASTA files with the genome reference sequences.
	-O <output dir>   
                (***) Specifies path to the directory of result files.
	-pass 
                Running STAR in the 2-pass alignment mode.
	-index
                Generating genome index of STAR with default settings.
	-SI_dir < genome index dir>  
                specifies path to the genome directory where genome indexes where generated.
	-I <path>  
                Specifies path to the directory where the eVIIS installation.
	-S <sample>   
                Name of sample
        -s <int>
                A parameter for featureCounts to indicate if strand-specific read counting should be performed. Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).

### Generating genome indexes ###
        Usage: python eVIIS.py -I /path2/eVIIS -index -SI_dir /path2/star_index -G /path2/refAnnoFile.gtf -GA /path2/refFastaFile.fasta
### Running STAR in the 2-pass mode, Evaluating virus expression level, Counting LS.Score and predicting TIME
        Usage: python eVIIS.py -I /path2/eVIIS -pass -SI_dir /path2/star_index -f1 seq1.fq.gz -f2 seq2.fq.gz -fq_dir /path2/seqs -G /path2/refAnnoFile.gtf -GA /path2/refFastaFile.fasta -O /path2/output -S sample -s strandSpecific

“””
Files Needed:
             refFastaFile.fasta
             refAnnoFile.gtf 
             seq1.fq.gz
             seq2.fq.gz
   
Output: 
       sample_mapped_reads.bam  
       sample_stringtie.gtf
       sample_virus.txt          (***)
       sample_featureCount.txt
       sample_count.txt
       sample_fpkm.txt
       sample_time.txt           (***)
“””

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####################################################### Predict the immune statuses (TIMEs) of multiple samples in a FPKM profile using eVIIS #####################################################
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
One step, Predicting immune statuses of multiple samples in a FPKM profile.
Usage: python predictTIME_FPKM_profile.py [OPTIONS]

The arguments of predictTIME_FPKM_profile.py are as followings:
	-I <path>  
                Specifies path to the directory where the eVIIS installation.
	-fp <path_to_FPKM_profile file>   
                Specifies the path to FPKM_profile.txt file needing to be analyzed.
        -OT_dir <dir of time files>  
                specifies path to the directory where the time files generated for all samples in FPKM_profile.txt file are.
        -OR <path_to_result_all file>
                specifies path to result_all.txt file which contains the TIME prediction result of all samples in FPKM_profile.txt file

### Running predictTIME_FPKM_profile.py to count LS.Scores and predict TIMEs for a FPKM_profile
        Usage: python predictTIME_FPKM_profile.py -I /path2/eVIIS -fp /path2/FPKM_profile.txt -OT_dir /path2/sample_time -OR /path2/result_all.txt

“””
Files Needed:
             FPKM_profile.txt 
             (Each row of FPKM_profile.txt is a gene, each column of FPKM_profile.txt is a sample. The first row is sample names, the first column is Gene-ids (like ENSG00000237851.1).)
   
Output: 
       sample1_time.txt, sample2_time.txt, …, samplen_time.txt
       result_all.txt           (***)
“””


The description of Output files generated by eVIIS:
---------------------------------------------------
1. The format of sample_virus.txt file
The top rows in this file are details about the transcripts expressed by virus(es) in corresponding sample, filtered from the StringTie result file.
*	-------------------- Virus expression level in this sample is: --------------------
*	virus: the name of virus detected in this sample.  
*	TPM: the expression level of corresponding virus in this sample. (In all transcripts expressed by the same virus, the one with the maximum TPM value was utilized to present the viral mRNA expression level.)

2. The format of sample_time.txt file
*	The first row shows the LF.Score of this sample.
*	The forth row shows the final TIME prediction result of this sample, as one of “Immune Exclusion”, “Immune Anergy”, “Immune Stimulation” 

3. The format of result_all.txt file
*	Sample: every sample name in FPKM_profile.txt file 
*	LF.Score: the LF.Score of corresponding sample calculated by eVIIS 
*	TIME: the immune status (TIME) of corresponding sample predicted by eVIIS

