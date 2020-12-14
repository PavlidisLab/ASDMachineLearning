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

* Guilt by association machine learning methods (GBA ML):

1) Princeton
    a) Paper: Krishnan et al., (2016) Genome-wide prediction and functional characterization of the genetic basis of autism spectrum disorder. Nat Neuro
    b) Data: 
        * training labels: Supplementary table 1
        * gene scores: Supplemtnary table 3
        * [http://asd.princeton.edu/] 
2) FRN
    a) Paper: Duda et al., (2018) Brain-specific functional relationship networks inform autism spectrum disorder gene prediction. Trans Psych.
    b) Data: [https://github.com/GuanLab/ASD_FRN]
    
3) DAMAGES
    a) Paper: Zhang and Shen (2017), A Cell Type-Specific Expression Signature Predicts Haploinsufficient Autism-Susceptibility Genes. Hum Mutat
    b) Data: 
        * gene scores + tranining labels: Supplemental table 4
    
4) RF_Lin
    a) Paper: Lin et al., (2018) A machine learning approach to predicting autism risk genes: Validation of known genes and discovery of new candidates
    b) Data: Supplementary table 2

5) PANDA
    a) Paper: Zhang et al., (2020), PANDA: Prioritization of autism-genes using network-based deep-learning approach. Genet Epidemiol
    b) Data: [https://github.com/MIB-Lab/PANDA]

* Hybrid GBA ML - genetics methods

1) forecASD
    a) Paper: Brueggeman et al., (2020) Forecasting risk gene discovery in autism with machine learning and genome-scale data. Sci Reports
    b) Data: [https://github.com/LeoBman/forecASD]
    
2) DAWN
    a) Paper: Liu et al., (2014) DAWN: a framework to identify autism genes and subnetworks using gene expression and genetics. Mol Autism
    b) Data: Supplemental table S4
    
* Genetic association methods (GA)

1) DeRubeis
    a) Paper: 
    b) Data: 
    
2) Sanders
    a) Paper: 
    b) Data: 
    
3) iHart
    a) Paper: 
    b) Data: 
    
4) Satterstrom
    a) Paper: 
    b) Data: 
    
5) Spark
    a) Paper: 
    b) Data: 
    
6) Iossifov
    a) Paper: 
    b) Data: 
    
* SFARI

1) SFARI
    b) Data: SFARIGene released 01-03-2020

* Constraint scores

1) ExAC pLI, mis_z
    a) Paper: Lek et al., (2016) Analysis of protein-coding genetic variation in 60,706 humans. Nature
    b) Data: [https://gnomad.broadinstitute.org/about]
    
2) gnomad pLI, oe_LoF
    a) Paper: Karczewski et al., (2020) The mutational constraint spectrum quantified from variation in 141,456 humans. Nature
    b) Data: [https://gnomad.broadinstitute.org/about]

* Generic scores

1) Number of publications
    b) Data: In-house calculation of number of PubMed publications per gene.
    
2) Multifunctionaltiy score
    b) Data: In-house calculation of number of functions per gene based on GO anntotation. 
    
3) Number of physical interaction partners
    b) Data: Calcualted from BioGrid version="3.5.169". 
