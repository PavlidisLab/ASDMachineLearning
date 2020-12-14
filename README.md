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

#### Guilt by association machine learning methods (GBA ML):

* *Princeton*
    * Paper: Krishnan et al., (2016) Genome-wide prediction and functional characterization of the genetic basis of autism spectrum disorder. Nat Neuro
    * Data: 
        * training labels: Supplemental table 1
        * gene scores: Supplemental table 3; probability
        * http://asd.princeton.edu/
* *FRN*
    * Paper: Duda et al., (2018) Brain-specific functional relationship networks inform autism spectrum disorder gene prediction. Trans Psych.
    * Data: https://github.com/GuanLab/ASD_FRN; supplemental table 1
        * gene score column (Supplemental table 1): Prediction Score 
    
* *DAMAGES*
    * Paper: Zhang and Shen (2017), A Cell Type-Specific Expression Signature Predicts Haploinsufficient Autism-Susceptibility Genes. Hum Mutat
    * Data: 
        * gene scores + tranining labels: Supplemental table 4
            * gene score column: Ensemble score
    
* *RF_Lin*
    * Paper: Lin et al., (2018) A machine learning approach to predicting autism risk genes: Validation of known genes and discovery of new candidates
    * Data: Supplemental table 2
        * gene score column: Mean

* *PANDA*
    * Paper: Zhang et al., (2020), PANDA: Prioritization of autism-genes using network-based deep-learning approach. Genet Epidemiol
    * Data: https://github.com/MIB-Lab/PANDA
        * gene score column: PredictionScore

####  Hybrid GBA ML - genetics methods

* *forecASD*
    * Paper: Brueggeman et al., (2020) Forecasting risk gene discovery in autism with machine learning and genome-scale data. Sci Reports
    * Data: https://github.com/LeoBman/forecASD
        * gene score column: forecASD
    
* *DAWN*
    * Paper: Liu et al., (2014) DAWN: a framework to identify autism genes and subnetworks using gene expression and genetics. Mol Autism
    * Data: Supplemental table 4
        * gene score column: min_lFDR
    
####  Genetic association methods (GA)

* *DeRubeis*
    * Paper: De Rubeis et al., (2014) Synaptic, transcriptional and chromatin genes disrupted in autism. Nature
    * Data: Supplemental table 2
        * gene score column: qvalue
    
* *Sanders*
    * Paper: Sanders et al., (2015) Insights into Autism Spectrum Disorder Genomic Architecture and Biology from 71 Risk Loci. Neuron. 
    * Data: Supplemental table 6
        * gene score column: tadaFdrAscSscExomeSscAgpSmallDel
    
* *iHart*
    * Paper: Ruzzo et al., (2019) Inherited and De Novo Genetic Risk for Autism Impacts Shared Networks. Cell
    * Data: Supplemental table 3
        * gene score column: qval
    
* *Satterstrom*
    * Paper: Satterstrom et al., (2020) Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism. Cell
    * Data: Supplemental tables 4, 5
        * gene score: 4; qval_dnccPTV
    
* *Spark*
    * Paper: Feliciano et al., (2019) Exome sequencing of 457 autism families recruited online provides evidence for autism risk genes. NPJ Genome Med
    * Data: Supplemental table 7
        *gene score column: QVal2
    
* *Iossifov*
    * Paper: Iossifov et al., (2015) Low load for disruptive mutations in autism genes and their biased transmission. PNAS
    * Data: Supplemental table 1, 2
        * gene score column: GR_PUBNOAUT_LGDs_post
    

#### SFARI

* SFARIGene released 01-03-2020

####  Constraint scores

* *ExAC pLI, mis_z*
    * Paper: Lek et al., (2016) Analysis of protein-coding genetic variation in 60,706 humans. Nature
    * Data: https://gnomad.broadinstitute.org/about
    
* *gnomad pLI, oe_LoF*
    * Paper: Karczewski et al., (2020) The mutational constraint spectrum quantified from variation in 141,456 humans. Nature
    * Data: https://gnomad.broadinstitute.org/about

#### Generic scores

* *Number of publications*
    * Data: In-house calculation of number of PubMed publications per gene.
    
* *Multifunctionaltiy score*
    * Data: In-house calculation of number of functions per gene based on GO anntotation. 
    
* *Number of physical interaction partners*
    * Data: Calcualted from BioGrid version="3.5.169". 
