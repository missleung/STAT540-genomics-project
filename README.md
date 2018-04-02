
**Studying the effect of DNA methylation on gene expression in human dorsolateral prefrontal cortex**
=====================================================================================================

This is the private repository of team **Gene Heroes** for STAT 540/ BIOF 540/ GSAT 540 term project. 


Members |  Department
---------|------------
Tiam Heydari | *Department of Biophysics*
Sina Jafarzadeh|	*Department of Bioinformatics*
Lisa Leung | *Department of Statisitcs*
Zohreh Sharafian|	*Department of Experimental Medicine*
Hiwot Tafessu |*Department of Statistics*


**Brief description of directories:**
---------------------------------------

**Step 1: Data Processing**:- Includes description of data, codes for data preprocessing in preparation for data analysis and image outputs from this step. 

**Step 2: Single Probe Analysis**:- This directory includes code and images on PCA analysis, eQTM analysis and correlation results of each probe and gene pair.

**Step 3: Multiple Probe Analysis**:- This directory will include the code and images on the exploration of multiple probe regression (full model using all probes that are significant), nonlinear regression, and LASSO variable selection, Stepwise (Forward and Backward) variable selection, and neural network method.

**Step 4: Biological Plausability**:- This directory includes the visualizations and explanation behind biological reasonings
and gene enrichment analysis. 

**shinyApp**:- Data visiualization app to supplement analysis outputs and visualize expression and methylation data for selected genes 

**zz_Deliverable**:- Project proposal, progress report and poster presentaion 


**Introduction**
----------------

DNA methylation is an epigenetic mechanism which plays a role in regulating tissue specific gene expression 1. DNA methylation mostly occurs when a cytosine base is next to a guanine base, forming a CpG site. Early studies showed that methylation of CpG sites prevent the expression of genes. Ho wever, recent studies revealed that methylation can be linked with both decreasing and increasing of gene expressions 2. In this study, we investigate the correlation between DNA methylation and gene expression in a quantitative manner (eQTM), i.e. evaluating the statistical significance of the effect of methylated CpG sites on gene expressions. In this project, we study the association of methylation and gene expression in human dorsolateral prefrontal cortex (DLPFC) 3.

**Goal**
----------------

There are a few studies quantizing the relationship of methylation and gene expression in a tissue specific manner4. However, this project has several advantages in terms of data type and analysis methods. 1) previous study investigate eQTM analysis in common tissue type including whole blood4, while we perform eQTM  of two data sets derived from the dorsolateral prefrontal cortex (DLPFC)5,6.2) we study the collective effect of multiple CpG sites on the expression of genes, while Bonder et al. only addressed the relationship of single CpG sites and genes. 3) we quantitatively assess the roles of chromatin states7, CpG site distance from gene and chromosome number on the regulatory relationship between methylation probes and genes. This provides us more biological intuition behind the regulatory processes in human brain.

**The Data**
----------------
We use a combination of two publicly available datasets, ROSMAP, namely as The Religious Orders Study (ROS)6 and The Memory and Aging Project (MAP)5. Subjects in ROSMAP were collected in 1994 and 1997 respectively by the same researchers. Both data sets are longitudinal cohort studies containing DNA methylation probes and gene expression data, and phenotypes of diseases related to aging. However, we do not use the phenotypes for the purposes of this project due to the confidentiality agreement. The gene expression data is RNA-sequencing data extracted from dorsolateral prefrontal cortex from 540 subjects using Illumina HiSeq with 101-bp paired-end reads. The DNA methylation data was generated using 450K Illumina array from 740 subjects. There are total of 468 samples individuals who have both the data. The data is QC’ed and corrected for common batch effect and co-founders 5,6.


**Summary of analysis and major results**
------------------------

.......


*References:*
---------------------------------------

........

