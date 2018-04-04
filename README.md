
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


**Task Division Brief description of directories:**
---------------------------------------

**[Step 1: Data Processing](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-1-Data%20Processing) (Sina)**:- Includes description of data, codes for data preprocessing in preparation for data analysis and image outputs from this step. 

**[Step 2: Single Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-2-Single%20Probe%20Analysis) (Lisa)**:- This directory includes code and images on PCA analysis, eQTM analysis and correlation results of each probe and gene pair. 

**[Step 3: Multiple Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-3-Multiple%20Probe%20Analysis) (Tiam and Sina)**:- This directory will include the code and images on the exploration of multiple probe regression (full model using all probes that are significant), nonlinear regression, and LASSO variable selection, Stepwise (Forward and Backward) variable selection, and neural network method. 

**[Step 4: Biological Plausability](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-4-BIological%20Plausibility) (Zohreh)**:- This directory includes the visualizations and explanation behind biological reasonings
and gene enrichment analysis. 

**[shinyApp](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/shinyApp) (Hiwot)**:- Data visiualization app to supplement analysis outputs and visualize expression and methylation data for selected genes 

**zz_Deliverable**:- [Project proposal](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Group_Project_Proposal_Gene_Heroes.md) (All), [progress report](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Progressreport.md) (All) and [poster presentation](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Poster_draft%20(Zohreh%20Sharafian's%20conflicted%20copy%202018-04-02%20(1))%20(Zohreh%20Sharafian's%20conflicted%20copy%202018-04-02)%20(tiam%20heydari's%20conflicted%20copy%202018-04-03).pdf) (Tiam, Zohreh, Hiwot, and Sina). 

**Introduction**
----------------

DNA methylation is an epigenetic mechanism which plays a role in regulating tissue specific gene expression 1. DNA methylation mostly occurs when a cytosine base is next to a guanine base, forming a CpG site. Early studies showed that methylation of CpG sites prevent the expression of genes. Ho wever, recent studies revealed that methylation can be linked with both decreasing and increasing of gene expressions 2. In this study, we investigate the correlation between DNA methylation and gene expression in a quantitative manner (eQTM), i.e. evaluating the statistical significance of the effect of methylated CpG sites on gene expressions. In this project, we study the association of methylation and gene expression in human dorsolateral prefrontal cortex (DLPFC) 3.

**Goal**
----------------

There are a few studies quantizing the relationship of methylation and gene expression in a tissue specific manner4. However, this project has several advantages in terms of data type and analysis methods. 1) previous study investigate eQTM analysis in common tissue type including whole blood4, while we perform eQTM  of two data sets derived from the dorsolateral prefrontal cortex (DLPFC)5,6.2) we study the collective effect of multiple CpG sites on the expression of genes, while Bonder et al. only addressed the relationship of single CpG sites and genes. 3) we quantitatively assess the roles of chromatin states7, CpG site distance from gene and chromosome number on the regulatory relationship between methylation probes and genes. This provides us more biological intuition behind the regulatory processes in human brain.

**The Data**
----------------
We use a combination of two publicly available datasets, [ROSMAP](https://www.synapse.org/#!Synapse:syn3219045), namely as The Religious Orders Study (ROS)6 and The Memory and Aging Project (MAP)5. Subjects in ROSMAP were collected in 1994 and 1997 respectively by the same researchers. Both data sets are longitudinal cohort studies containing DNA methylation probes and gene expression data, and phenotypes of diseases related to aging. However, we do not use the phenotypes for the purposes of this project due to the confidentiality agreement. The gene expression data is RNA-sequencing data extracted from dorsolateral prefrontal cortex from 540 subjects using Illumina HiSeq with 101-bp paired-end reads. The DNA methylation data was generated using 450K Illumina array from 740 subjects. There are total of 468 samples individuals who have both the data. The data is QC’ed and corrected for common batch effect and co-founders 5,6.


**Summary of analysis and major results**
------------------------

***[Step 1: Data Processing](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1.md)***(https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1.md): ROSMAP raw data is preprocessed to prepare for the single probe and multiple probe analyses. We started with a gene expression data set, DNA methylation data set and its corresponding CpG sites. For this particular assignment, we are only interested in the probes within a 1Mb distance around the gene of interest. 

In order to demonstrate our multiple probe analysis, we have to handle with the issue of multicollinearity between probes.  Hence, we decide to pick one probe out of a cluster of probes that are strongly correlated with each other (ie. probes that have strong linkage disequilibrium). We used A-clust algorithm to solve this issue. Based on the resulting plots of [number of non-singletone clusters](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-1.png) and [average intra-cluster correlation](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-2.png), we decide that correlation threshold of 0.8 and a base-pair distance of 1000 are the ideal parameters to separate the probes. In return, we are left with three data sets for Step 2 and 3: [probes_genes_distance, probes_subjects and subjects_genes](https://drive.google.com/file/d/1Ieze9KwSy5UL9c6Vt5uFF9ta4U9kvkoe/).

***[Step 2: Single Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2.md)***: Using the [data processed in Step 1](https://drive.google.com/file/d/1Ieze9KwSy5UL9c6Vt5uFF9ta4U9kvkoe/) , the PCA is used to adjust for any hidden confounders that could exist in our data sets (in both gene expression and DNA methylation). To choose the number of PCs adjusted for both gene expression data and DNA methylation probe data, we plotted the [heatmap](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/prop_of_genes_sig_on_50_genes_V4_fdr.png) on the proportion of genes that have at least one significant probe with FDR < 0.1 based on number of PCs adjusted for both data sets. Note that the heatmap is based on 50 genes. For exploration, we also plotted the percentage of cumulative variance that additional PCs are adjusted for in both [gene expression data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/cum_var_explained_50_V4_gdata.png) and [DNA methylation data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/cum_var_explained_50_V4_probe.png). The optimal number of PCs adjusted will be the highest number of significant values based on the heatmap. 

As a [result](https://drive.google.com/open?id=1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY), we concluded that the optimal DNA methylation probe data set should be adjusted for five PCs and gene expression data set is adjusted for 30 PCs. To observe [the effects of adjusting for number of PCs in the gene expression data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/gene_C21orf56_gdata.png), we took a look at a gene that show the most number of significant probes. Similarly, [the effects of adjusting for number of PCs in DNA methylation data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/gene_C21orf56_probes.png) are explored in the same gene.  


***[Step 3: Multiple Probe Analysis on Statistical Models](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3.md)***: Based on the single probe analysis, we would like to explore what happens when we modelled linear regression with more than one probe for each gene epxression. Specifically, we ran a single linear regression (using the probe with the smallest FDR value) and compared with the multiple linear regression using all significant probes (FDR < 0.05). In addition, we also ran Forward and Backward stepwise variable selection method for each gene, LASSO variable selection method and 

***[Step 3: Multiple Probe Analysis on Neural Network]()***

***[Step 4: Biological Plausability](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4.md)***: Based on the single probe analysis in Step 1, we looked at the distances of probe and genes which are correlated at FDR < 0.05 in a [plot](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4_files/figure-markdown_github/unnamed-chunk-7-1.png). We observe that most of these pairs are very close by the distance, with very few pairs that are far from each other. Based on another [plot](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4_files/figure-markdown_github/unnamed-chunk-8-1.png), we also see that the average FDR increases when distances from probe and gene pair widen.jb

**Discussion and Conclusion**
----------------

*References:*
---------------------------------------

........

