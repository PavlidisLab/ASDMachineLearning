# ASDMachineLearning

## Intro

Discovering genes involved in complex human genetic disorders is a major challenge. Many have suggested that machine learning (ML) algorithms using gene networks can be used to supplement traditional genetic association-based approaches to predict or prioritize disease genes. However, questions have been raised about the utility of ML methods for this type of task due to biases within the data, and poor real-world performance. Using autism spectrum disorder (ASD) as a test case, we sought to investigate the question: Can machine learning aid in the discovery of disease genes? We collected thirteen published ASD gene prioritization studies and evaluated their performance using known and novel high-confidence ASD genes. We also investigated their biases towards generic gene annotations, like number of association publications. We found that ML methods which do not incorporate genetics information have limited utility for prioritization of ASD risk genes. These studies perform at a comparable level to generic measures of likelihood for the involvement of genes in any condition, and do not out-perform genetic association studies. Future efforts to discover disease genes should be focused on developing and validating statistical models for genetic association, specifically for association between rare variants and disease, rather than developing complex machine learning methods using complex heterogeneous biological data with unknown reliability.

Preprint: [Can machine learning aid in identifying disease genes?](https://doi.org/10.1101/2020.11.26.394676)

Publication: [“Guilt by association” is not competitive with genetic association for identifying autism risk genes](https://www.nature.com/articles/s41598-021-95321-y)

## Code

All scripts are written in the R programming language: 

*allClassifiers*:
* Performance Evaluation: allClassifiers_evaluate.R 
* Confidence Interval Plots: allClassifiers_CI_plots.R
* Correlation Plots: allClassifiers_correlate_plot.R

*forecASD retest*:
* forecASD Model Modification: forecASD_03_ensemble_model.adaptations.R


## Session Info

```{r, echo=FALSE}
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] zoo_1.8-8           DescTools_0.99.38   ggpubr_0.4.0        tidyr_1.1.2         purrr_0.3.4         ggplot2_3.3.2      
[7] dplyr_1.0.2         randomForest_4.6-14

loaded via a namespace (and not attached):
 [1] textshape_1.7.1   tidyselect_1.1.0  reshape2_1.4.4    lattice_0.20-38   haven_2.3.1       carData_3.0-4     expm_0.999-5     
 [8] colorspace_1.4-1  vctrs_0.3.4       generics_0.0.2    rlang_0.4.8       e1071_1.7-4       pillar_1.4.6      foreign_0.8-71   
[15] glue_1.4.2        withr_2.3.0       readxl_1.3.1      rootSolve_1.8.2.1 lifecycle_0.2.0   plyr_1.8.6        stringr_1.4.0    
[22] munsell_0.5.0     ggsignif_0.6.0    gtable_0.3.0      cellranger_1.1.0  zip_2.1.1         mvtnorm_1.1-1     rio_0.5.16       
[29] forcats_0.5.0     lmom_2.8          curl_4.3          class_7.3-15      broom_0.7.1       Rcpp_1.0.5        backports_1.1.10 
[36] scales_1.1.1      abind_1.4-5       gld_2.6.2         Exact_2.1         hms_0.5.3         stringi_1.5.3     openxlsx_4.2.2   
[43] rstatix_0.6.0     grid_3.6.0        tools_3.6.0       magrittr_1.5      tibble_3.0.4      crayon_1.3.4      car_3.0-10       
[50] pkgconfig_2.0.3   Matrix_1.2-17     ellipsis_0.3.1    MASS_7.3-51.4     data.table_1.13.0 rstudioapi_0.11   R6_2.4.1         
[57] boot_1.3-22       compiler_3.6.0  
```

