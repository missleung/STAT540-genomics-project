Gene Heros: Step1- Preprocessing and Data Preparation
================
Sina Jafarzadeh
January 13, 2018

``` r
library(ggplot2)
library(pheatmap)
library(magrittr)
library(reshape2)
```

***1.1 Data Introduction***

Our data consists of the gene expression and methylation measurements for 508 and 702 subjects, respectively. We organized the data into the following files:

-   gene\_expression.csv: a 420103\*702 matrix showing the values for ~420 K CpG sites across 702 samples.
-   methylation.csv: a 508\*13484 matrix showing the values for ~13 K genes across 508 samples.
-   genes\_name.csv: a 13484 elements vector containing the gene symbol and the ensemble name for all genes.
-   cpg\_sites\_name: a 420103 elements vector containing the names for all cpg sites according to the Illumina Infinium HumanMethylation450K BeadChip naming convention.
-   subjects\_in\_gene\_expression\_data.csv: the cols in gene expression matrix presenting the name of subjects that we collected gene expression data for.
-   subjects\_in\_methylation\_data.csv: the rows in gene expression matrix presenting the name of subjects that we collected gene expression data for.

Apparently, it imposes an unnecessary huge number of hypothesis, testing the correlation between each single CpG site and gene pair. So, We solely analyze the probes locating in a 1MB window around the gene of interest using the information encoded in probe\_gene\_distance\_matrix\_sparse.csv.

***1.2 Multicollinearity of Methylation Probes***

It is known that methylation probes can be in strong Linkage Disequilibrium (LD) with their neighbouring probes. So, we need to reduce the number of probes to a set of unrelated ones. It helps us to decrease the computational cost of the model. Moreover, It is necessary to have uncorrelated probes in the second phase of the project where we want to provide a linear predictor of multiple probes for the gene of interest. There are a few tools developed for this purpose. In this project, we use A-clust algorithm provided by Harvard University. This algorithm first merges each pair of nearby probes and all of the probes located between them if the correlation between that pair is more than corr\_dist\_thresh. Nearby pairs are pairs of probes with a base pair distance less than bp\_dist\_thresh. In the next step, it merges each two neighbour probes if their correlation distance is less than corr\_dist\_thresh. We find the best values for these parameters by creating train and validation datasets. We choose the parameters based on train data and then see if they work on another validation data. We choose CpG sites on chromosomes 1-18 as train data and CpG sites in chromosome 19-22 as validation. We choose the parameters that decrease the number of probes (i.e. create as much cluster as possible) and yield the highest intra-cluster correlation between probes. Below, we plot the results for train set chromosomes:



    report_matrix_clusts_number_non_singles_train <- read.table("a-clustering/report_matrix_clusts_number_non_singles_train.csv", header=FALSE, sep=",")
    report_matrix_intra_cluster_corr_non_singles_train <- read.table("a-clustering/report_matrix_corrs_non_singles_train.csv", header=FALSE, sep=",")

    row_col_combinaions = expand.grid(corr_threshold=c("0.9", "0.8", "0.6", "0.2"), bp_dist_threshold=c("500","1000","2000","4000"))
    report_matrix_clusts_number_non_singles_train_melt = report_matrix_clusts_number_non_singles_train %>% melt()
    report_matrix_intra_cluster_corr_non_singles_train_melt = report_matrix_intra_cluster_corr_non_singles_train %>%  melt()

    heatmap_data = cbind(row_col_combinaions, number_of_clusters=report_matrix_clusts_number_non_singles_train_melt$value)
    heatmap_data %>% ggplot(aes(corr_threshold,bp_dist_threshold,fill = number_of_clusters)) + geom_tile() + scale_fill_gradient("Number of non-singleton clusters",low = 'blue',high = 'red') +xlab('Correlation threshold')+ylab('BP distance threshold') 

    heatmap_data = cbind(row_col_combinaions, intra_cluster_correlation=report_matrix_intra_cluster_corr_non_singles_train_melt$value)
    heatmap_data %>% ggplot(aes(corr_threshold,bp_dist_threshold,fill = intra_cluster_correlation)) + geom_tile() + scale_fill_gradient("Average intra-cluser correlation",low = 'blue',high = 'red') +xlab('Correlation threshold')+ylab('BP distance threshold') 

The plots suggest that corr\_threshold = 0.8 and bp\_dist\_threshold = 1000 are the best parameters in terms of the intra-cluster correlation/number of clusters trade-off. The authors of A-clust paper also suggested this parameter in their experiment on human methylation data, supporting the chosen parameters in our problem. Now, we investigate the parameters on validation set. Below, you see the similar plots for validation data:

``` r
report_matrix_clusts_number_non_singles_test <- read.table("../Results/a-clust/report_matrix_clusts_number_non_singles_test.csv", header=FALSE, sep=",")
report_matrix_intra_cluster_corr_non_singles_test <- read.table("../Results/a-clust/report_matrix_corrs_non_singles_test.csv", header=FALSE, sep=",")

row_col_combinaions = expand.grid(corr_threshold=c("0.9", "0.8", "0.6", "0.2"), bp_dist_threshold=c("500","1000","2000","4000"))
report_matrix_clusts_number_non_singles_test_melt = report_matrix_clusts_number_non_singles_test %>% melt()
```

    ## No id variables; using all as measure variables

``` r
report_matrix_intra_cluster_corr_non_singles_test_melt = report_matrix_intra_cluster_corr_non_singles_test %>%  melt()
```

    ## No id variables; using all as measure variables

``` r
heatmap_data = cbind(row_col_combinaions, number_of_clusters=report_matrix_clusts_number_non_singles_test_melt$value)
heatmap_data %>% ggplot(aes(corr_threshold,bp_dist_threshold,fill = number_of_clusters)) + geom_tile() + scale_fill_gradient("Number of non-singleton clusters",low = 'blue',high = 'red') +xlab('Correlation threshold')+ylab('BP distance threshold') 
```

![](Step1_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
heatmap_data = cbind(row_col_combinaions, intra_cluster_correlation=report_matrix_intra_cluster_corr_non_singles_test_melt$value)
heatmap_data %>% ggplot(aes(corr_threshold,bp_dist_threshold,fill = intra_cluster_correlation)) + geom_tile() + scale_fill_gradient("Average intra-cluser correlation",low = 'blue',high = 'red') +xlab('Correlation threshold')+ylab('BP distance threshold') 
```

![](Step1_files/figure-markdown_github/unnamed-chunk-1-2.png)

It verifies the consistency of parameters fitness in train and validation data. The final results are in aclust\_sampled folder. We would need the following files:

-   methylationSNMnorm\_usr\_prb\_matrix\_sampled\_corr\_dist\_thresh\_point2\_bp\_dist\_thresh\_1000\_train.csv

-   methylationSNMnorm\_usr\_prb\_matrix\_sampled\_corr\_dist\_thresh\_point2\_bp\_dist\_thresh\_1000\_test.csv

***1.3 Data Convertion***

In this part, we convert the described .csv files to R language data format to make it easier to handle. We also combine the related .csv files to organize the data better. In order to lower the computational cost, we limit the analysis to chromosomes 19-22. Below is the list of final data structure saves together as rosmap\_postdocV1.RData file:

-   probes\_genes\_distance: a sparse matrix created by Matrix library showing the distance of each cpg probe to each gene in its neighbourhood.
-   probes\_subjects: a data frame showing the Z-scored value for each cpg site across all subjects.
-   subjects\_genes: a data frame showing the expression value of each gene across all subjects.

We can access the list of probes, genes, and subjects using the rownames() and colnames() functions.