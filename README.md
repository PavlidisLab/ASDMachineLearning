# ASDMachineLearning

## Intro

Discovering genes involved in complex human genetic disorders is a major challenge. Many have suggested that machine learning (ML) algorithms using gene networks can be used to supplement traditional genetic association-based approaches to predict or prioritize disease genes. However, questions have been raised about the utility of ML methods for this type of task due to biases within the data, and poor real-world performance. Using autism spectrum disorder (ASD) as a test case, we sought to investigate the question: Can machine learning aid in the discovery of disease genes? We collected thirteen published ASD gene prioritization studies and evaluated their performance using known and novel high-confidence ASD genes. We also investigated their biases towards generic gene annotations, like number of association publications. We found that ML methods which do not incorporate genetics information have limited utility for prioritization of ASD risk genes. These studies perform at a comparable level to generic measures of likelihood for the involvement of genes in any condition, and do not out-perform genetic association studies. Future efforts to discover disease genes should be focused on developing and validating statistical models for genetic association, specifically for association between rare variants and disease, rather than developing complex machine learning methods using complex heterogeneous biological data with unknown reliability.

Preprint: Can machine learning aid in identifying disease genes? The case of autism spectrum disorder
URL: [https://doi.org/10.1101/2020.11.26.394676](https://doi.org/10.1101/2020.11.26.394676)


## Code

All scripts are written in the R programming language: 

* Performance Evalaution: allClassifiers_evaluate.R 
* Confidence Internval Plots: allClassifiers_CI_plots.R
* Correlation Plots: allClassifiers_correlate_plot.R

* forecASD Model Modification: 


## Packages

## Session Info

## Source Data

