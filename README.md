
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

**[Step 1: Data Processing](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-1-Data%20Processing)** - Includes description of data, codes for data preprocessing in preparation for data analysis and image outputs from this step. __Sina contributed in both the research and coding in Step 1.__ 

**[Step 2: Single Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-2-Single%20Probe%20Analysis)** - This directory includes code and images on PCA analysis, eQTM analysis and correlation results of each probe and gene pair. __Hiwot contributed in researching the topic and Lisa contributed in the coding.__

**[Step 3: Multiple Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-3-Multiple%20Probe%20Analysis)** - This directory will include the code and images on the exploration of multiple probe regression (full model using all probes that are significant), nonlinear regression, and LASSO variable selection, Stepwise (Forward and Backward) variable selection, and neural network method. __Tiam contributed researching and coding using regression modelling and variable selection methods for Step 3 (full model, nonlinear regression, Lasso and Stepwise methods), and Sina contributed in the researching and coding in neural network method.__

**[Step 4: Biological Plausability](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-4-BIological%20Plausibility)** - This directory includes the visualizations and explanation behind biological reasonings and gene enrichment analysis. __Sina contributed in the research and idea of gene enrichment analysis and Zohreh contributed in the coding and written aspect of the gene enrichment analysis.__

**[shinyApp](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/shinyApp)** - Data visiualization app to supplement analysis outputs and visualize expression and methylation data for selected genes. __Hiwot coded the shinyApp program.__ 

**[zz_Deliverable](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/zz_Deliverable)**:

- [Project proposal](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Group_Project_Proposal_Gene_Heroes.md) - __Lisa prepared the proposal, Hiwot and Tiam contributed to the ideas, then Zohreh and Sina finalized the proposal.__
- [Progress report](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Progressreport.md) - __Lisa prepared the report, Hiwot and Tiam contributed to the ideas, then Zohreh and Sina finalized the report.__
- [Poster presentation](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Poster_draft%20(Zohreh%20Sharafian's%20conflicted%20copy%202018-04-02%20(1))%20(Zohreh%20Sharafian's%20conflicted%20copy%202018-04-02)%20(tiam%20heydari's%20conflicted%20copy%202018-04-03).pdf) - __Zohreh and Sina prepared and finalized the poster, then Hiwot and Tiam provided organized all the graphs plots.__

**Repo Preparation** - Organizing and preparing the entire repository at presentable level. __Hiwot and Lisa organized the Repo. Zohreh and Sina wrote Introduction, Background and Data Description. Lisa prepared the summary analyses but are finalized by everyone accordingly to the steps they have contributed. Zohreh wrote Bibliography.__

**Introduction**
----------------

DNA methylation is an epigenetic mechanism which plays a role in regulating tissue specific gene expression 1. DNA methylation mostly occurs when a cytosine base is next to a guanine base, forming a CpG site. Early studies showed that methylation of CpG sites prevent the expression of genes. Ho wever, recent studies revealed that methylation can be linked with both decreasing and increasing of gene expressions 2. In this study, we investigate the correlation between DNA methylation and gene expression in a quantitative manner (eQTM), i.e. evaluating the statistical significance of the effect of methylated CpG sites on gene expressions. In this project, we study the association of methylation and gene expression in human dorsolateral prefrontal cortex (DLPFC) 3.

**Goal**
----------------

There are a few studies quantizing the relationship of methylation and gene expression in a tissue specific manner4. However, this project has several advantages in terms of data type and analysis methods. 1) previous study investigate eQTM analysis in common tissue type including whole blood4, while we perform eQTM  of two data sets derived from the dorsolateral prefrontal cortex (DLPFC)5,6.2) we study the collective effect of multiple CpG sites on the expression of genes, while Bonder et al. only addressed the relationship of single CpG sites and genes. 3) we quantitatively assess the roles of chromatin states7, CpG site distance from gene and chromosome number on the regulatory relationship between methylation probes and genes. This provides us more biological intuition behind the regulatory processes in human brain.


![](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/figone.png)


**The Data**
----------------
We use a combination of two publicly available datasets, [ROSMAP](https://www.synapse.org/#!Synapse:syn3219045), namely as The Religious Orders Study (ROS)6 and The Memory and Aging Project (MAP)5. Subjects in ROSMAP were collected in 1994 and 1997 respectively by the same researchers. Both data sets are longitudinal cohort studies containing DNA methylation probes and gene expression data, and phenotypes of diseases related to aging. However, we do not use the phenotypes for the purposes of this project due to the confidentiality agreement. The gene expression data is RNA-sequencing data extracted from dorsolateral prefrontal cortex from 540 subjects using Illumina HiSeq with 101-bp paired-end reads. The DNA methylation data was generated using 450K Illumina array from 740 subjects. There are total of 468 samples individuals who have both the data. The data is QCâed and corrected for common batch effect and co-founders 5,6.


**Summary of analysis and major results**
------------------------

**[Step 1: Data Processing](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1.md)**: ROSMAP raw data is preprocessed to prepare for the single probe and multiple probe analyses. We started with a gene expression data set, DNA methylation data set and its corresponding CpG sites. For this particular assignment, we are only interested in the probes within a 1Mb distance around the gene of interest. 

In order to demonstrate our multiple probe analysis, we have to handle with the issue of multicollinearity between probes.  Hence, we decide to pick one probe out of a cluster of probes that are strongly correlated with each other (ie. probes that have strong linkage disequilibrium). We used A-clust algorithm to solve this issue. Based on the resulting plots of [number of non-singletone clusters](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-1.png) and [average intra-cluster correlation](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-2.png), we decide that correlation threshold of 0.8 and a base-pair distance of 1000 are the ideal parameters to separate the probes. In return, we are left with three data sets for Step 2 and 3: [probes_genes_distance, probes_subjects and subjects_genes](https://drive.google.com/file/d/1Ieze9KwSy5UL9c6Vt5uFF9ta4U9kvkoe/).

**[Step 2: Single Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2.md)**: Using the [data processed in Step 1](https://drive.google.com/file/d/1Ieze9KwSy5UL9c6Vt5uFF9ta4U9kvkoe/) , the PCA is used to adjust for any hidden confounders that could exist in our data sets (in both gene expression and DNA methylation). To choose the number of PCs adjusted for both gene expression data and DNA methylation probe data, we plotted the [heatmap](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/prop_of_genes_sig_on_50_genes_V4_fdr.png) on the proportion of genes that have at least one significant probe with FDR < 0.1 based on number of PCs adjusted for both data sets. Note that the heatmap is based on 50 genes. For exploration, we also plotted the percentage of cumulative variance that additional PCs are adjusted for in both [gene expression data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/cum_var_explained_50_V4_gdata.png) and [DNA methylation data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/cum_var_explained_50_V4_probe.png). The optimal number of PCs adjusted will be the highest number of significant values based on the heatmap. 

As a [result](https://drive.google.com/open?id=1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY), we concluded that the optimal DNA methylation probe data set should be adjusted for five PCs and gene expression data set is adjusted for 30 PCs. To observe [the effects of adjusting for number of PCs in the gene expression data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/gene_C21orf56_gdata.png), we took a look at a gene that show the most number of significant probes. Similarly, [the effects of adjusting for number of PCs in DNA methylation data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/gene_C21orf56_probes.png) are explored in the same gene.  


**[Step 3: Multiple Probe Analysis on Statistical Models](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3.md)**: Based on the single probe analysis, we would like to explore what happens when we modeled linear regression with more than one probe for each gene expression. To achieve this, we first divided our data set to 80% train and 20% test. Then we ran a single linear regression (using the probe with the smallest FDR value) and compared with the multiple linear regression using all significant probes (FDR < 0.05). In addition, we also ran Forward and Backward stepwise variable selection method for each gene, LASSO variable selection method. to assess the effect of higher degree fitting on our data, we also checked the higher degree polynomial fits on the single probes.

The [results shows that the full probes set](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3_files/figure-markdown_github/unnamed-chunk-15-1.png) is the best model for more than 50% of the genes. interestingly the second best used the most significant probes transformed by a high degree function.
To deal with the danger of overfitting we also [plot the ratio of test to train error](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3_files/figure-markdown_github/unnamed-chunk-16-6.png) and we found that it is symmetrically distributed. around 1 which indicate that we didn't overfit the data.

**[Step 3: Multiple Probe Analysis on Neural Network]()**

**[Step 4: Biological Plausability](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4.md)**: Based on the single probe analysis in Step 1, we assess [the spatial distribution of significant methylation probes](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4-spatial_distribution_of_methylated_CpG_sites.md) according to their distance from genes at FDR < 0.05. The [histogram](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4_files/figure-markdown_github/unnamed-chunk-6-1.png) shows that the majority of methylated CpG sites are concentrated near the transcription start sites (TSS) of genes. The [bar graph](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4_files/figure-markdown_github/unnamed-chunk-9-1.png)shows that the methylated CpG sites with significant P-values are located close to TSS. This result is consistent with previous data which showed that CpG sites are relatively enriched around TSS of genes (Numata et al., 2012; Numata, Ye, Herman, & Lipska, 2014; Saxonov et al., 2006).
We also analyze the distribution of significant methylation probes in different chromatin states using ChromHMM (Ernst & Kellis, 2012). The bars show the frequency of total number of positive and negative regulating CpG sites for different chromatin states. It shows that TSS enriches with the highest number of significant methylation probes compared to other chromatin states.
Moreover, we see that the number of negative regulating methylation probes in TSS region is higher compared to the positive ones. Biological experiments also verify that most of the methylated CpG sites near the TSS and gene promotor
repress the expression value of the gene resulting in a down-regulation effect.
We also perform [“Gene Set enrichment analysis”](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/geneset_enrichment_analysis.md) to find the relationship between the number of significant probes and the complexity of gene regulation using erminR package (Gillis et al., 2010). The analysis shows that there is not
any significant relationship between the numbers of significant methylated CpG sites corresponding to genes and the multifunctionality.


**Discussion and Conclusion**
----------------

*References:*
---------------------------------------

........

