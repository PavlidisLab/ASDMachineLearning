###
#EVALUATE all classifiers
###


###
#libraries
###
library(dplyr)
library(DescTools)
library(purrr)


######

unzip("./Data/allClassifiers_forecASD_202012.zip")
allClassifiers_ROC = read.csv("./Data/allClassifiers_forecASD_202012.csv", stringsAsFactors = FALSE) %>%
    dplyr::select(-X) 



###
#rates / AUCs / ROCs / PRs
#evaluativeRates:
#prep function for getAUCs: 
#1)filters allClassifiers by evaluation set to remove training labels for point estimation and confidence intervals
#2) produces an ordered dataframe used in getACUS

#getACUS:
#used to get AURROC, P20R, P43R
#output is a single row dataframe per score


#score = name of score being evaluated = output will be ordered by this score: highest 1=ASD to 0=no ASD
#genesInPhenotype = gene set being evalauted = new.SFARI-HC or TADA_novel
#positiveLabels = training genes for score = i.e. forecASD_posLabs or "" for GA studies
#SFARIset = specify SFARIHC, new.SFARIHC to make sure the genes are removed if TADA_novel
#sampleSet = empty excpet for calcualting CIs

###

evaluativeRates = function(score, genesInPhenotype, positiveLabels, SFARIset, sampleSet){
    #score = "forecASD_score"
    #genesInPhenotype = "TADA_novel"
    #positiveLabels = "forecASD_posLabs"
    #SFARIset="new.SFARIHC"
    #sampleSet=""
    method.name = score
    
    #only used for bootstrapping: filters allClassifiers by the randomly sampled gene set made up of randomly sampled evaluation genes (newSFARIHC or TADA_novel) and the rest of the genes in the classifier
    if(length(sampleSet)>1) {
        geneList = sampleSet %>% unique
        name = geneList[1]
        geneList = geneList[2:length(geneList)]
        
        rankedGeneList = allClassifiers_ROC %>% dplyr::filter(primary.gene.symbol %in% geneList)
        genesInPhenotype_edit = genesInPhenotype
        
    } else{
        rankedGeneList = allClassifiers_ROC
        genesInPhenotype_edit = genesInPhenotype
        name = genesInPhenotype_edit
        
    }
        
    
    #no positive lables = used as the control evaluation = used for SFARIHC/new.SFARIHC eval
    if(positiveLabels==""){
        rankedGeneList = rankedGeneList[,c("primary.gene.symbol",score, genesInPhenotype_edit)] %>%
            dplyr::rename_at(., genesInPhenotype_edit, list(~gsub(genesInPhenotype_edit, "labels", .))) 
        
    }else {
        #used for novel evaluation
        if(grepl("TADA_novel", genesInPhenotype)){
            if(positiveLabels=="SFARIHC"|positiveLabels=="new.SFARIHC"){
                #used for GA/generic methods to remove "positive training genes" i.e. known ASD genes
                rankedGeneList = rankedGeneList[,c("primary.gene.symbol",score, genesInPhenotype_edit, positiveLabels)] %>%
                    dplyr::rename_at(., positiveLabels, list(~gsub(positiveLabels, "positiveLabels", .))) %>%
                    dplyr::filter(positiveLabels!=1)%>%
                    dplyr::rename_at(., genesInPhenotype_edit, list(~gsub(genesInPhenotype_edit, "labels", .))) 
                
            } else {
                if(SFARIset=="new.SFARIHC"){
                    #used for GBA ML methods to remove 1) new.SFARIHC genes AND 2) positive training genes
                    rankedGeneList = rankedGeneList[,c("primary.gene.symbol",score, genesInPhenotype_edit, positiveLabels, "new.SFARIHC")] %>%
                        dplyr::rename_at(., positiveLabels, list(~gsub(positiveLabels, "positiveLabels", .))) %>%
                        dplyr::filter(positiveLabels!=1)%>%
                        dplyr::filter(new.SFARIHC!=1)%>%
                        dplyr::rename_at(., genesInPhenotype_edit, list(~gsub(genesInPhenotype_edit, "labels", .)))
                    
                } else {
                    #used for GBA ML methods to remove 1) SFARIHC genes AND 2) positive training genes
                    rankedGeneList = rankedGeneList[,c("primary.gene.symbol",score, genesInPhenotype_edit, positiveLabels, "SFARIHC")] %>%
                        dplyr::rename_at(., positiveLabels, list(~gsub(positiveLabels, "positiveLabels", .))) %>%
                        dplyr::filter(positiveLabels!=1)%>%
                        dplyr::filter(SFARIHC!=1)%>%
                        dplyr::rename_at(., genesInPhenotype_edit, list(~gsub(genesInPhenotype_edit, "labels", .))) 
                    
                }
               
            }
        }
        
    }
    size = dim(rankedGeneList %>% dplyr::filter(labels==1))[1]
    
    p = sum(rankedGeneList$labels==1) #total number of true positives
    n = sum(rankedGeneList$labels==0) # total number of true negatives
    
    #oe lof is in the opposite scale, so have to rank the score the other way
    if(score=="oe_lof_upper_score"){
        rankedGeneList = rankedGeneList %>%
            dplyr::rename_at(.,vars(contains("_score")), list(~gsub(".*_", "", .))) %>% 
            dplyr::arrange(score)
        
    } else {
        #rank the gene score so 1= ASD is at the top and 0=not ASD is at the bottom
        rankedGeneList = rankedGeneList %>%
            dplyr::rename_at(.,vars(contains("_score")), list(~gsub(".*_", "", .))) %>% 
            dplyr::arrange(-score)
    }
    
    
    #calculate rates: AT EACH GENE SCORE THREHOLD = TP is the number of actual positives (labels=1=evaluation gene sets) and FP being everything that is not an actual positive
    #calculate TPR, FPR and Precision for EACH SCORE THRESHOLD as number of true predicted positives over total p, number of false predicted positives over total negatives, and predicted true positives over predicted true positives and predicted false positives
    rates= rankedGeneList %>%
        group_by(score) %>%
        dplyr::mutate(TP = sum(labels),
                      FP = sum(!labels)) %>%
        ungroup()%>%
        dplyr::select(score, TP, FP) %>%
        unique %>% 
        dplyr::mutate(TPR = cumsum(TP)/(p),
                      FPR = cumsum(FP)/(n),
                      Precision = cumsum(TP)/(cumsum(TP)+cumsum(FP))) %>%
        dplyr::mutate(TPC = cumsum(TP)) 
    
    #TruePosRate = cumsum(rankedGeneList$labels)/sum(rankedGeneList$labels)
    #FalsePosRate = cumsum(!rankedGeneList$labels)/sum(!rankedGeneList$labels)
    #Precision = cumsum(rankedGeneList$labels)/(cumsum(rankedGeneList$labels)+cumsum(!rankedGeneList$labels))
    
    dat = left_join(rates, rankedGeneList, by="score") %>% dplyr::mutate(size=size) %>% dplyr::mutate(name = name) %>% dplyr::mutate(method.name = method.name)
    
    dat = as.data.frame(dat)
    
    return(dat)
    
    
}

getAUCs = function(score, genesInPhenotype, positiveLabels, SFARIset, sampleSet){
    #score = "forecASD_score"
    #genesInPhenotype = "TADA_top_common"
    #positiveLabels = "SFARIHC"
    #controlGeneList = ""
    #sampleSet = ""
    
    #call evaluative rates
    rates = evaluativeRates(score, genesInPhenotype, positiveLabels, SFARIset, sampleSet)
    
    #get unique values for score thresholds
    ratesEdit = rates %>% dplyr::select(score, TPC, TPR, FPR, Precision) %>% unique
    
    #calcualte AUROC
    AUC = AUC((ratesEdit$FPR), (ratesEdit$TPR), "trapezoid")
    #calculate precision at 20% recall of positive labels
    Pr.20Recallnum = round(max(ratesEdit$TPC*.20))
    Pr.20Recall =  ratesEdit %>%
        dplyr::filter(TPC==Pr.20Recallnum) %>%
        dplyr::select(Precision) %>%
        dplyr::summarize(Pr.20Recall = mean(Precision)*100) %>%
        dplyr::mutate(Pr.20Recall = round(Pr.20Recall, 4))
    #calculate precision at 20% recall of positive labels if there is not an exact match between TPC and recall number i.e. for TIES and it skips over the number
    if(is.na(Pr.20Recall$Pr.20Recall)){
        Pr.20Recall = ratesEdit %>%
            dplyr::mutate(TPC_edit = abs(TPC-Pr.20Recallnum)) %>% 
            dplyr::filter(TPC_edit == min(TPC_edit)) %>%
            dplyr::select(Precision) %>%
            dplyr::summarize(Pr.20Recall = mean(Precision)*100) %>% 
            dplyr::mutate(Pr.20Recall =paste0(round(Pr.20Recall,4), "*"))
    }
    Pr.43Recallnum = round(max(ratesEdit$TPC*.43))
    Pr.43Recall = ratesEdit %>%
        dplyr::filter(TPC==Pr.43Recallnum) %>%
        dplyr::select(Precision) %>%
        dplyr::summarize(Pr.43Recall = mean(Precision)*100) %>%
        dplyr::mutate(Pr.43Recall = round(Pr.43Recall, 4))
    
    if(is.na(Pr.43Recall$Pr.43Recall)){
        Pr.43Recall = ratesEdit %>%
            dplyr::mutate(TPC_edit = abs(TPC-Pr.43Recallnum)) %>% 
            dplyr::filter(TPC_edit == min(TPC_edit)) %>%
            dplyr::select(Precision) %>%
            dplyr::summarize(Pr.43Recall = mean(Precision)*100) %>% 
            dplyr::mutate(Pr.43Recall =paste0(round(Pr.43Recall,4), "*"))
    }
    
    
    dat = data.frame(method=score, 
                     geneList = genesInPhenotype,
                     AUC=(AUC), 
                     Pr.20Recall=as.character(Pr.20Recall$Pr.20Recall), 
                     labelNum20=as.character(Pr.20Recallnum), 
                     Pr.43Recall=as.character(Pr.43Recall$Pr.43Recall), 
                     labelNum43=Pr.43Recallnum, 
                     geneListSize = unique(rates$size),
                     stringsAsFactors = FALSE)
    
    return(dat)

}


#test2 = getAUCs("forecASD_score", "new.SFARIHC", "","","")

#t3=Sys.time()
#test4 = getAUCs("forecASD_score", "TADA_novel", "forecASD_posLabs","new.SFARIHC","")
#t4=Sys.time()


#############
#bootstrap

#bootstrap for getAUCS for 95% CIS:

#score = name of score being evaluated = output will be ordered by this score: highest 1=ASD to 0=no ASD
#genesInPhenotype = gene set being evaluated = new.SFARI-HC or TADA_novel
#positiveLabels = training genes for score = i.e. forecASD_posLabs or "" for GA studies
#SFARIset = specify SFARIHC, new.SFARIHC to make sure the genes are removed if TADA_novel

#1) calls evaluativeRates to get ordered dataframe
#2) randomly samples all genes 
#2.a) pool1 = lables i.e. evalautive gene set SFARIHC or TADA_novel
#2.b) pool2 = non-labeles i.e. all the rest of the genes in allClassifiers
#2.c) randomly sample 2500 times from each pool with replacement, and combine the pool1 and pool2 into one sampleSet
#3) sampleSet is fed into getAUCs 
#3.1) evaluativeRates will: 
    #1) keep the UNIQUE genes from each sampleSet
    #2) filter allClassifiers by the RANDOMLY SAMPLED geneSet made up of pool 1 and pool 2 i.e. so a randomly sampled set of both the evaluation gene set, and the rest of the genes in the classifier is left over
    #3) and note which of the geneInPhenotype (evaluation gene set, new.SFARIHC, or TADA_novel genes) are still left in the random sample to use in getAUC
    #4) and filter out positive training labels if necessary
    #5) rank the filtered allClassifiers by the selected score
#4) getAUCs will:
    #1) evaluate the randomly sampled classifier using the genesInPhenotype (evluation gene set, new.SFHARIHC or TADA_novel
    #2) will output the AUROC, P20R etc.
    #3) will repeat 2500 times, and output the 95% confidence intervals

#############

evaluativeAUCRatesBootstrap = function(score, genesInPhenotype, positiveLabels, SFARIset){
    #score = "forecASD_score"
    #genesInPhenotype = "TADA_novel"
    #positiveLabels = "forecASD_posLabs"
    #SFARIset= "new.SFARIHC"
    
    #realRates = evaluativeRates(score, genesInPhenotype, positiveLabels, "", "")
    #realAUCs = getAUCs(score, genesInPhenotype, positiveLabels, "", "")
    
    #realRates_edit = realRates %>% dplyr::select(primary.gene.symbol, labels)
    
    #set seeed for replication
    set.seed(50)
    
    #point estimate
    realRates = evaluativeRates(score, genesInPhenotype, positiveLabels, SFARIset,"")
    realAUCs = getAUCs(score, genesInPhenotype, positiveLabels, SFARIset, "")
    
    #spilt by labels (genes =1 = evalation gene set; genes=0=everything else), to randomly sample by each to create a new gene set to evaluate
    pool1 = realRates %>% dplyr::select(primary.gene.symbol, labels) %>% dplyr::filter(labels==1)
    pool2 = realRates %>% dplyr::select(primary.gene.symbol, labels) %>% dplyr::filter(labels==0)
    
    #sample 2500 times from each pool
    sam1 = replicate(n=2500, sample(pool1$primary.gene.symbol, size = nrow(pool1), replace = TRUE))
    colnames(sam1)=paste0("sampleSet", 1:ncol(sam1))
    sam2 = replicate(n=2500, sample(pool2$primary.gene.symbol, size = nrow(pool2), replace = TRUE))
    colnames(sam2)=paste0("sampleSet", 1:ncol(sam2))
    #and combine so the new gene set will be a mash-up of positive genes and the rest of the genes in the classifier
    sampleSet = rbind(sam1, sam2)
    
    #masterset with all 2500 sample (control) geneSets
    sampleSet = split(sampleSet, rep(1:ncol(sampleSet), each = nrow(sampleSet)))
    #names(controlGeneList) = paste0(genesInPhenotype, "control_", names(controlGeneList))
    
    #feed all 2500 sample gene sets through getAUCS
    #1) goes through evaluativeRates to get a NEW ORDERED LIST that has ONLY THE RANDOMLY SAMPLED GENES
    #2) goes through getAUCs to get the vals ONLY USING THE POSITIVE LABELED GENES FROM THE RANDOM SAMPLE THAT ARE IN TEH ORDERED LIST
    arglist = list(score = rep(score, 2500),
                   genesInPhenotype = rep(genesInPhenotype,2500),
                   positiveLabels = rep(positiveLabels, 2500),
                   SFARIset=SFARIset,
                   sampleSet = sampleSet)
    
    t1 = Sys.time()
    boot_sam = pmap_dfr(arglist, getAUCs)
    t2= Sys.time()
    
    boot_sam = boot_sam %>%
        dplyr::mutate(Pr.20RecallE = gsub("\\*", "", Pr.20Recall) %>% as.numeric,
                      Pr.43RecallE = gsub("\\*", "", Pr.43Recall)  %>% as.numeric)
    
    #95% CI by quantiles
    quantileAUC = quantile(boot_sam$AUC, probs = c(0.025, 0.975))
    quantilePR20 = quantile(boot_sam$Pr.20RecallE, probs = c(0.025, 0.975))
    quantilePr43 = quantile(boot_sam$Pr.43RecallE, probs = c(0.025, 0.975))
    
    dat = data.frame(score = score,
                     geneList = genesInPhenotype,
                     realAUC = round(realAUCs$AUC,5),
                     realP20R = realAUCs$Pr.20Recall,
                     realP43R = realAUCs$Pr.43Recall,
                     AUC_CI = toString(c(round(quantileAUC[["2.5%"]],5), round(quantileAUC[["97.5%"]],5))),
                     PR20_CI =toString(c(round(quantilePR20[["2.5%"]],5), round(quantilePR20[["97.5%"]],5))),
                     PR43_CI = toString(c(round(quantilePr43[["2.5%"]],5), round(quantilePr43[["97.5%"]],5))))
    return(dat)
    
}


####
#4 evaluations:
#1) all scores on the new SFARIHC (2020) genes (cnrl expt)
#scores = all 15 = 7 GBA ML; 5 GA; 3 constraint
#genesInPhenotype = all new.SFARIHC
#positiveLabels = "" = empty b/c all built using these gene
#SFARIset = "" = empty = only for TADA_novel

####

t1 = Sys.time()
arglist = list(score = c("ASDprinceton_score", "ASD_frn_score", "DAMAGES_score","RF_Lin_score", "forecASD_score", "DAWN_score","PANDA_score",
                         "DeRubeis_score", "Sanders_score", "iHart_score","Satterstrom_score", "Iossifov_score",
                         "exac_pLI_score", "gnomad_pLI_score", "oe_lof_upper_score"),
               genesInPhenotype = rep("new.SFARIHC", 15),
               positiveLabels = rep("",15),
               SFARIset=rep("",15))
new.SFARIHCbootstrapaucs = pmap_dfr(arglist, evaluativeAUCRatesBootstrap)
write.csv(new.SFARIHCbootstrapaucs, "./Results/SFARIHC_bootstrap_202012.csv")
t2= Sys.time()

####
#2 evaluations:
#2) all scores on the TADA_novel and new.SFARIHC genes = test expt
#scores = all 15 = 7 GBA ML; 5 GA; 3 constraint
#genesInPhenotype = all TADA_novel
#positiveLabels = 
    #GBA ML and HYBRID training genes = "ASDprinceton_posLabs", "ASDfrn_posLabs", "DAMAGES_posLabs", "ASDfrn_posLabs", "forecASD_posLabs", "new.SFARIHC" (DAWN DOESN'T HAVE TRAINING GENES),"PANDA_posLabs",
    #GA and generic = don't have training genes; remove SFARIHC genes; "new.SFARIHC", "new.SFARIHC", "new.SFARIHC","new.SFARIHC", "new.SFARIHC","new.SFARIHC","new.SFARIHC", "new.SFARIHC" 
#SFARIset = new.SFARIHC = removes these genes

####


t1 = Sys.time()
arglist = list(score = c("ASDprinceton_score", "ASD_frn_score", "DAMAGES_score","RF_Lin_score", "forecASD_score", "DAWN_score","PANDA_score",
                         "DeRubeis_score", "Sanders_score", "iHart_score","Satterstrom_score", "Iossifov_score",
                         "exac_pLI_score", "gnomad_pLI_score", "oe_lof_upper_score"),
               genesInPhenotype = rep("TADA_novel", 15),
               positiveLabels = c("ASDprinceton_posLabs", "ASDfrn_posLabs", "DAMAGES_posLabs", "ASDfrn_posLabs", "forecASD_posLabs", "new.SFARIHC","PANDA_posLabs",
                                  "new.SFARIHC", "new.SFARIHC", "new.SFARIHC","new.SFARIHC", "new.SFARIHC","new.SFARIHC","new.SFARIHC", "new.SFARIHC"),
               SFARIset=rep("new.SFARIHC",15))
new.NOVELbootstrapaucs = pmap_dfr(arglist, evaluativeAUCRatesBootstrap)
write.csv(new.NOVELbootstrapaucs, "./Results/TADAnovel_bootstrap_202012.csv")

t2= Sys.time()




####
#4 evaluations:
#3) modified forecASD scores on the new SFARIHC (2020) genes (cnrl expt)
#scores = all 6 forecASD models = redo_forec_score, noClass_forec_score,noClassPPI_forec_score,noClassPPIBS_forec_score, PPIonly_forec_score,BrainSpanOnly_forec_score
#genesInPhenotype = all new.SFARIHC
#positiveLabels = "" = empty b/c all built using these gene
#SFARIset = "" = empty = only for TADA_novel
####
arglist = list(score = c("redo_forec_score", "noClass_forec_score","noClassPPI_forec_score",
                         "noClassPPIBS_forec_score", "PPIonly_forec_score","BrainSpanOnly_forec_score"),
               genesInPhenotype = rep("new.SFARIHC", 6),
               positiveLabels = rep("",6),
               SFARIset= rep("",6))
new.SFARIHCbootstrapaucs = pmap_dfr(arglist, evaluativeAUCRatesBootstrap)
write.csv(new.SFARIHCbootstrapaucs, "./Results/forecASD_SFARIHC_bootstrap_202012.csv")


####
#4 evaluations:
#4) modified forecASD scores on the TADA_novel genes, and removed positive training lables/SFARIHC genes (test expt)
#scores = all 6 forecASD models = redo_forec_score, noClass_forec_score,noClassPPI_forec_score,noClassPPIBS_forec_score, PPIonly_forec_score,BrainSpanOnly_forec_score
#genesInPhenotype = all TADA_novel
#positiveLabels = 
#forecASD training genes =  "forecASD_posLabs",
#SFARIset = new.SFARIHC = removes these genes

####


arglist =list(score = c("redo_forec_score", "noClass_forec_score","noClassPPI_forec_score",
                        "noClassPPIBS_forec_score", "PPIonly_forec_score","BrainSpanOnly_forec_score"),
              genesInPhenotype = rep("TADA_novel", 6),
              positiveLabels = rep("forecASD_posLabs",6),
              SFARIset= rep("new.SFARIHC",6))
new.TADAnovelbootstrapaucs = pmap_dfr(arglist, evaluativeAUCRatesBootstrap)
write.csv(new.TADAnovelbootstrapaucs, "./Results/forecASD_TADAnovel_bootstrap_202012.csv")


####

