# ASDMachineLearning
## Scripts

All scripts are written in the R programming language: 

*allClassifiers*:
* Performance Evalaution: allClassifiers_evaluate.R 
    * Calcualting 2500 stratified bootstraped samples to obtain 95% confidence intervals of AUROC and Precision at 20% recall of ASD gene sets
    
* Confidence Internval Plots: allClassifiers_CI_plots.R
    * AUC and PR graphs of GBA ML, Hybrid and GA methods on ASD gene sets
    * 95% confidence intervals of AUROC and Precision and 20% recall of ASD gene sets by  GBA ML, Hybrid and GA methods
    
* Correlation Plots: allClassifiers_correlate_plot.R
    * Spearman correlation heatmap of all GBA ML, Hybrid, GA, Constraint, and Generic gene annotations

*forecASD retest*:
* forecASD Model Modification: forecASD_03_ensemble_model.adaptations.R
    * Modified from https://github.com/LeoBman/forecASD
    * Rerun to obtain final randomForest models using different feature sets to ascertain performance on recovering ASD gene sets, and feature importance
