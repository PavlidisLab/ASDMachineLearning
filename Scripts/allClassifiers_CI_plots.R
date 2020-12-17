###
#PLOT CIS
###

###
#libraries
###
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(ggpubr)

#####
unzip("./Data/allClassifiers_forecASD_202012.zip")
allClassifiers_ROC = read.csv("./Data/allClassifiers_forecASD_202012.csv", stringsAsFactors = FALSE) %>%
    dplyr::select(-X) 


####
#CI PLOTS and Tables
####
new.SFARIHC_CI = read.csv( "./Results/SFARIHC_bootstrap_202012.csv", stringsAsFactors = FALSE) %>%  dplyr::select(-X) %>%dplyr::mutate(SFARIset="")
new.novel_CI = read.csv( "./Results/TADAnovel_bootstrap_202012.csv", stringsAsFactors = FALSE)%>%  dplyr::select(-X) %>%dplyr::mutate(SFARIset="new.SFARIHC")

forecASD.new.SFARIHC_CI = read.csv( "./Results/forecASD_SFARIHC_bootstrap_202012.csv", stringsAsFactors = FALSE) %>%  dplyr::select(-X) %>%dplyr::mutate(SFARIset="")
forecASD.new.novel_CI = read.csv( "./Results/forecASD_TADAnovel_bootstrap_202012.csv", stringsAsFactors = FALSE)%>%  dplyr::select(-X) %>%dplyr::mutate(SFARIset="new.SFARIHC")

allCIs = rbind( new.SFARIHC_CI ,new.novel_CI,forecASD.new.SFARIHC_CI,forecASD.new.novel_CI) %>%
    dplyr::rename(., "method.name"="score")



####
#CI PLOTS and Tables

#type = GBA_ML, GA, Generic
#score = name of score being evaluated = output will be ordered by this score: highest 1=ASD to 0=no ASD
#genesInPhenotype = gene set being evalauted = new.SFARI-HC or TADA_novel
#positiveLabels = training genes for score = i.e. forecASD_posLabs or "" for GA studies
#SFARIset = specify SFARIHC, new.SFARIHC to make sure the genes are removed if TADA_novel
#sampleSet = empty excpet for calcualting CIs

####
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


plotROC_PR = function(type,score, genesInPhenotype, positiveLabels, SFARIset, sampleSet ){
    
    #type="GBA_ML"
    #score = c("ASDprinceton_score", "ASD_frn_score","DAMAGES_score", "RF_Lin_score","PANDA_score", "forecASD_score", "DAWN_score")
    #genesInPhenotype = rep("TADA_novel", 7)
    #positiveLabels = c("ASDprinceton_posLabs", "ASDfrn_posLabs","DAMAGES_posLabs", "ASDfrn_posLabs","PANDA_posLabs","forecASD_posLabs","new.SFARIHC")
    #SFARIset = rep("new.SFARIHC",7)
    #sampleSet= rep("",7)
    
    arglist = list(score, genesInPhenotype, positiveLabels, SFARIset, sampleSet)
    
    e_rates = pmap_dfr(arglist, evaluativeRates)
    
    
    if(type=="GBA_ML"){
        
        
        levels1 = c("Princeton (PRI)", "FRN", "DAMAGES (DAM)","RF_Lin (RFL)", "PANDA (PAN)", "forecASD (FOR)", "DAWN (DAW)")
        levels2 = c("PRI", "FRN", "DAM","RFL","PAN", "FOR", "DAW")
        
        colours = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
        
        #title1 = paste(paste("","        GBA ML", sep="\t      "), "A)", sep="\n")
        #title1 = paste("\t\t\tGBA ML\nA)") 
        title2 = "B)"
        title3 = "C)"
        title4 = "D)"
        
        
        if(unique(genesInPhenotype)=="TADA_novel"){
            title1 = paste(paste("    ", "    GBA ML on novel-HC", sep="\t\t"), "A)", sep="\n")
            ylabA = ("ROC")
            ylabB = ("AUROC 95% CI")
            ylabC = ("Precision Recall ")
            ylabD = ("P20R 95% CI")
        } else {
            title1 = paste(paste("    ", "    GBA ML on SFARI-HC", sep="\t\t"), "A)", sep="\n")
            ylabA = ("ROC")
            ylabB = ("AUROC 95% CI")
            ylabC = ("Precision Recall")
            ylabD = ("P20R 95% CI")
            
            add.Precision = e_rates %>% 
                dplyr::group_by(method.name) %>% 
                dplyr::filter(score==max(score)) %>% 
                dplyr::select(method.name,score, Precision) %>% 
                unique
            
            dat = data.frame(score=rep(NA,7),
                             TP=rep(NA,7),
                             FP = rep(NA,7),
                             TPR=rep(0,7),
                             FPR=rep(NA,7),
                             Precision = add.Precision$Precision,
                             TPC = rep(NA,7),
                             method.name=add.Precision$method.name)
            e_rates = e_rates %>% dplyr::select(score, TP, FP, TPR, FPR, Precision,TPC, method.name) %>% rbind(., dat)
            
        }
        
        
    } else if(type=="Genetics"){
        
        levels1 = c("DeRubeis (DER)", "Sanders (SAN)", "iHart", "Satterstrom (SAT)", "Iossifov (IOS)")
        levels2 = c("DER", "SAN", "iHart", "SAT", "IOS")
        
        colours = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3" ,"#FF7F00")
        if(unique(genesInPhenotype)=="TADA_novel"){
            title1 = paste("Genetic Association on novel-HC\nI)") 
        } else {
            title1 = paste("Genetic Association on SFARI-HC\nI)") 
        }
        #title1 = paste(paste("","     Genetic Association", sep="\t"), "I)", sep="\n")
        #title1 = paste("\t\tGenetic Association\nI)") 
        title2 ="J)"
        title3 = "K)"
        title4 = "L)"
        
        ylabA = ""
        ylabB = ""
        ylabC =""
        ylabD = ""
        
        
        
        
        add.Precision = e_rates %>% 
            dplyr::group_by(method.name) %>% 
            dplyr::filter(score==max(score)) %>% 
            dplyr::select(method.name,score, Precision) %>% 
            unique
        
        
        dat = data.frame(score=rep(NA,5),
                         TP=rep(NA,5),
                         FP = rep(NA,5),
                         TPR=rep(0,5),
                         FPR=rep(NA,5),
                         Precision = add.Precision$Precision,
                         TPC = rep(NA,5),
                         method.name=add.Precision$method.name)
        e_rates = e_rates %>% dplyr::select(score, TP, FP, TPR, FPR, Precision,TPC, method.name) %>% rbind(., dat)
        
        
        
    } else if(type=="Generic") {
        
        levels1 = c("ExAC_pLI", "gnomAD_pLI", "oe_LoF")
        levels2 = c("ExAC_pLI", "gnomAD_pLI", "oe_LoF")
        
        colours = c("#386CB0","#F0027F","#BF5B17")
        
        #title1 = paste(paste("","     Generic", sep="\t\t  "), "E)", sep="\n")
        #title1 = paste("\t\tGeneric\nE)") 
        if(unique(genesInPhenotype)=="TADA_novel"){
            title1 = paste(paste("    ", "    Generic on novel-HC", sep="\t\t"), "E)", sep="\n")
        } else {
            title1 = paste(paste("    ", "    Generic on SFARI-HC", sep="\t\t"), "E)", sep="\n")
        }
        
        title2 = "F)"
        title3 ="G)"
        title4 = "H)"
        
        ylabA = ""
        ylabB = ""
        ylabC =""
        ylabD = ""
        
        add.Precision = e_rates %>% 
            dplyr::group_by(method.name) %>% 
            dplyr::filter((grepl("pLI", method.name)&score==max(score))|grepl("lof", method.name)&score==min(score)) %>% 
            dplyr::select(method.name,score, Precision) %>% 
            unique
        
        dat = data.frame(score=rep(NA,3),
                         TP=rep(NA,3),
                         FP = rep(NA,3),
                         TPR=rep(0,3),
                         FPR=rep(NA,3),
                         Precision = add.Precision$Precision,
                         TPC = rep(NA,3),
                         method.name=add.Precision$method.name)
        e_rates =  e_rates %>% dplyr::select(score, TP, FP, TPR, FPR, Precision,TPC, method.name) %>% rbind(., dat)
        
        
        
    } else {
        levels1 = c("Redo", "NoClass (NoC)", "NoClassPPI (NoCP)", "NoClassPPIBS (NoCPB)", "PPIOnly (PPI)", "BrainSpanOnly (BS)")
        levels2 = c("Redo", "NoC", "NoCP", "NoCPB", "PPI", "BS")
        
        #colours = c("#3288BD", "#5E4FA2", "#1B7837","#C51B7D" ,"#8C510A","#35978F")
        #colours = c("#8DD3C7", "#FFFFB3" ,"#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462")
        #colours = c("#7EBDB2", "#c2c232", "#AAA7C2", "#CF6A5F", "#638BA3", "#C28A4C")
        #colours = c("#1F78B4", "#33A02C","#E31A1C","#FF7F00" ,"#6A3D9A","#B15928")
        #colours = c("#A6CEE3","#B2DF8A","#FB9A99" ,"#FDBF6F","#CAB2D6","#FFFF99")
        colours = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")
        
        col = "forecASD Version"
        x = "forecASD Version"
        
        
        if(unique(genesInPhenotype)=="TADA_novel"){
            title1 = paste(paste("","forecASD on novel-HC genes", sep="\t\t\t"), "A)", sep="\n")
            title2 = "B)"
            title3 = "C)"
            title4 = "D)"
            
            ylabA = ("ROC")
            ylabB = ("AUROC 95% CI")
            ylabC = ("Precision Recall")
            ylabD = ("P20R 95% CI")
        } else {
            title1 = paste(paste("","forecASD on SFARI-HC genes", sep="\t\t\t"), "E)", sep="\n")
            title2 = "F)"
            title3 = "G)"
            title4 = "H)"
            
            ylabA = ("ROC")
            ylabB = ("AUROC 95% CI")
            ylabC = ("Precision Recall")
            ylabD = ("P20R 95% CI")
            
            add.Precision = e_rates %>% 
                dplyr::group_by(method.name) %>% 
                dplyr::filter(score==max(score)) %>% 
                dplyr::select(method.name,score, Precision) %>% 
                unique
            
            dat = data.frame(score=rep(NA,6),
                             TP=rep(NA,6),
                             FP = rep(NA,6),
                             TPR=rep(0,6),
                             FPR=rep(NA,6),
                             Precision = add.Precision$Precision,
                             TPC = rep(NA,6),
                             method.name=add.Precision$method.name)
            e_rates = e_rates %>% dplyr::select(score, TP, FP, TPR, FPR, Precision,TPC, method.name) %>% rbind(., dat)
            
        }
    }
    
    
    
    
    g1=e_rates %>%
        dplyr::select(score, TPR, FPR, Precision, method.name)%>%
        dplyr::filter(score %in% score) %>%
        dplyr::mutate(method.name = gsub("_score", "", method.name),
                      method.name = ifelse(method.name=="exac_pLI", "ExAC_pLI",
                                           ifelse(method.name=="gnomad_pLI", "gnomAD_pLI",
                                                  ifelse(method.name=="oe_lof_upper", "oe_LoF", method.name)))) %>%
        dplyr::mutate(method.name = ifelse(method.name=="ASDprinceton", "Princeton (PRI)",
                                           ifelse(method.name=="ASD_frn", "FRN",
                                                  ifelse(method.name=="DAMAGES", "DAMAGES (DAM)",
                                                         ifelse(method.name=="DAWN", "DAWN (DAW)",
                                                                ifelse(method.name=="RF_Lin", "RF_Lin (RFL)",
                                                                       ifelse(method.name=="forecASD", "forecASD (FOR)",
                                                                              ifelse(method.name=="PANDA", "PANDA (PAN)", method.name))))))))%>%
        dplyr::mutate(method.name=ifelse(method.name=="DeRubeis", "DeRubeis (DER)",
                                         ifelse(method.name=="Sanders", "Sanders (SAN)",
                                                ifelse(method.name=="Satterstrom", "Satterstrom (SAT)",
                                                       ifelse(method.name=="Iossifov", "Iossifov (IOS)", method.name)))))%>%
        dplyr::mutate(method.name = ifelse(method.name == "redo_forec", "Redo",
                                           ifelse(method.name=="noClass_forec", "NoClass (NoC)",
                                                  ifelse(method.name=="noClassPPI_forec", "NoClassPPI (NoCP)",
                                                         ifelse(method.name=="noClassPPIBS_forec", "NoClassPPIBS (NoCPB)",
                                                                ifelse(method.name=="PPIonly_forec", "PPIOnly (PPI)",
                                                                       ifelse(method.name=="BrainSpanOnly_forec", "BrainSpanOnly (BS)", method.name))))))) %>%
        dplyr::mutate(method.name= factor(method.name, levels=levels1))%>%
        #dplyr::mutate(method.name= factor(method.name))%>%
        ggplot(aes(x=TPR, y=Precision, col=method.name))+
        geom_line(size=1) +
        #geom_abline(slope=1,intercept=0, size=1, col="grey")+
        xlim(0,1)+
        ylim(0,1)+
        scale_colour_manual(labels = levels1,
                            values=colours)+
        scale_fill_discrete(name = "", labels = levels1)+
        labs(#title = paste0("C) ", title1, "Precision Recall for", title2),
            title = title3,
            col = "Paper")+
        ylab(paste(ylabC, "Precision", sep="\n"))+
        theme_bw()+
        guides(col=guide_legend(nrow=3))+
        theme(legend.title = element_blank())
    #theme(plot.title = element_text(size=34, family="Helvetica"))+
    #theme(axis.text = element_text(size=26, family="Helvetica"),
    #      axis.title = element_text(size=28, family="Helvetica"),
    #      legend.text = element_text(size=24, family="Helvetica"),
    #      legend.title = element_text(size=2, family="Helvetica"))
    
    g2=e_rates %>%
        dplyr::select(score, TPR, FPR, Precision, method.name)%>%
        rbind(dat = data.frame(score = c(NA,NA),
                               TPR=c(0,1),
                               FPR=c(0,1),
                               Precision =c(0,1),
                               method.name="")) %>%
        dplyr::filter(score %in% score) %>%
        dplyr::mutate(method.name = gsub("_score", "", method.name),
                      method.name = ifelse(method.name=="exac_pLI", "ExAC_pLI",
                                           ifelse(method.name=="gnomad_pLI", "gnomAD_pLI",
                                                  ifelse(method.name=="oe_lof_upper", "oe_LoF", method.name)))) %>%
        dplyr::mutate(method.name = ifelse(method.name=="ASDprinceton", "Princeton (PRI)",
                                           ifelse(method.name=="ASD_frn", "FRN",
                                                  ifelse(method.name=="DAMAGES", "DAMAGES (DAM)",
                                                         ifelse(method.name=="DAWN", "DAWN (DAW)",
                                                                ifelse(method.name=="RF_Lin", "RF_Lin (RFL)",
                                                                       ifelse(method.name=="forecASD", "forecASD (FOR)", 
                                                                              ifelse(method.name=="PANDA", "PANDA (PAN)",method.name))))))))%>%
        dplyr::mutate(method.name=ifelse(method.name=="DeRubeis", "DeRubeis (DER)",
                                         ifelse(method.name=="Sanders", "Sanders (SAN)",
                                                ifelse(method.name=="Satterstrom", "Satterstrom (SAT)",
                                                       ifelse(method.name=="Iossifov", "Iossifov (IOS)", method.name)))))%>%
        dplyr::mutate(method.name = ifelse(method.name == "redo_forec", "Redo",
                                           ifelse(method.name=="noClass_forec", "NoClass (NoC)",
                                                  ifelse(method.name=="noClassPPI_forec", "NoClassPPI (NoCP)",
                                                         ifelse(method.name=="noClassPPIBS_forec", "NoClassPPIBS (NoCPB)",
                                                                ifelse(method.name=="PPIonly_forec", "PPIOnly (PPI)",
                                                                       ifelse(method.name=="BrainSpanOnly_forec", "BrainSpanOnly (BS)", method.name))))))) %>%
        
        dplyr::mutate(method.name= factor(method.name, levels=c(levels1, "")))%>%
        ggplot(aes(x=FPR, y=TPR, col=method.name))+
        geom_line(size=1)+
        #geom_abline(slope=1,intercept=0, size=1, col="grey")+
        xlim(0,1)+
        ylim(0,1)+
        scale_color_manual(breaks = c(levels1),
                           values=c(colours, "#696969"))+
        
        labs(#title = paste0("A) ",title1, "ROC for", title2),
            title = title1,
            col = "Paper")+
        ylab(paste(ylabA, "TPR", sep="\n"))+
        theme_bw()+
        guides(col=guide_legend(nrow=3))+
        theme(legend.title = element_blank())
    #theme(plot.title = element_text(size=34, family="Helvetica"))+
    #theme(axis.text = element_text(size=26, family="Helvetica"),
    #      axis.title = element_text(size=28, family="Helvetica"),
    #      legend.text = element_text(size=24, family="Helvetica"),
    #      legend.title = element_text(size=2, family="Helvetica"),
    #      legend.position="none")
    
    
    #####
    #CI: GBA_ML
    #####
    
    
    g3=allCIs %>%
        dplyr::filter(method.name %in% score) %>%
        dplyr::filter(geneList %in% unique(genesInPhenotype)) %>%
        dplyr::rename(., "SFARIsetFilteredOut"="SFARIset")%>%
        dplyr::filter(SFARIsetFilteredOut %in% unique(SFARIset))%>%
        dplyr::mutate(method.name = gsub("_score", "", method.name),
                      method.name = ifelse(method.name=="exac_pLI", "ExAC_pLI",
                                           ifelse(method.name=="gnomad_pLI", "gnomAD_pLI",
                                                  ifelse(method.name=="oe_lof_upper", "oe_LoF", method.name))))%>%
        dplyr::mutate(method.name = ifelse(method.name=="ASDprinceton", "PRI",
                                           ifelse(method.name=="ASD_frn", "FRN",
                                                  ifelse(method.name=="DAMAGES", "DAM",
                                                         ifelse(method.name=="DAWN", "DAW",
                                                                ifelse(method.name=="RF_Lin", "RFL",
                                                                       ifelse(method.name=="forecASD", "FOR",
                                                                              ifelse(method.name=="PANDA", "PAN", method.name))))))))%>%
        dplyr::mutate(method.name=ifelse(method.name=="DeRubeis", "DER",
                                         ifelse(method.name=="Sanders", "SAN",
                                                ifelse(method.name=="Satterstrom", "SAT",
                                                       ifelse(method.name=="Iossifov", "IOS", method.name))))) %>%
        dplyr::mutate(method.name = ifelse(method.name == "redo_forec", "Redo",
                                           ifelse(method.name=="noClass_forec", "NoC",
                                                  ifelse(method.name=="noClassPPI_forec","NoCP",
                                                         ifelse(method.name=="noClassPPIBS_forec", "NoCPB",
                                                                ifelse(method.name=="PPIonly_forec", "PPI",
                                                                       ifelse(method.name=="BrainSpanOnly_forec", "BS", method.name))))))) %>%
        dplyr::mutate(method.name = factor(method.name, levels = levels2))%>%
        
        separate(., AUC_CI, c("p5AUC", "p95AUC"), ", ") %>% 
        separate(., PR20_CI, c("p5PR", "p95PR"), ", ") %>%
        dplyr::mutate(p5AUC = as.numeric(p5AUC), p95AUC=as.numeric(p95AUC))%>%
        dplyr::mutate(p5PR = as.numeric(p5PR),
                      p95PR=as.numeric(p95PR), 
                      realP20R = as.numeric(gsub("\\*$", "", realP20R))) %>%
        #dplyr::mutate(colour = brewer.pal(7, "Dark2"))%>%
        dplyr::select(method.name, realP20R , p5PR, p95PR) %>%
        unique%>%
        ggplot(aes(method.name, realP20R))+
        geom_point(aes(colour=method.name), size=1)+
        #facet_wrap(~category,scales="free_x", nrow=1)+
        ylim(0,100)+
        geom_errorbar(aes(ymin=p5PR, ymax=p95PR, width=.2),
                      position = position_dodge(0.05), size=1,colour=colours) +
        scale_colour_manual(values = colours) +
        xlab("")+
        labs(#title = paste0("D) ", title1, "P20R 95% Confidence Intervals"),
            title = title4,
            col = "Paper",
            x="Paper")+
        ylab(paste(ylabD, "P20R (%)", sep="\n"))+
        theme_bw()+
        guides(col=guide_legend(nrow=3))+
        theme(legend.text = element_text(size=2))+
        theme(legend.title = element_blank())
    #theme(plot.title = element_text(size=34, family="Helvetica"))+
    #theme(axis.text = element_text(size=26, family="Helvetica"),
    #      axis.title = element_text(size=28, family="Helvetica"),
    #      legend.text = element_text(size=24, family="Helvetica"),
    #      legend.title = element_text(size=2, family="Helvetica"),
    #      legend.position="none")
    
    g4=allCIs %>%
        dplyr::filter(method.name %in% score) %>%
        dplyr::filter(geneList %in% unique(genesInPhenotype)) %>%
        dplyr::rename(., "SFARIsetFilteredOut"="SFARIset")%>%
        dplyr::filter(SFARIsetFilteredOut %in% unique(SFARIset))%>%
        dplyr::mutate(method.name = gsub("_score", "", method.name),
                      method.name = ifelse(method.name=="exac_pLI", "ExAC_pLI",
                                           ifelse(method.name=="gnomad_pLI", "gnomAD_pLI",
                                                  ifelse(method.name=="oe_lof_upper", "oe_LoF", method.name)))) %>%
        dplyr::mutate(method.name = ifelse(method.name=="ASDprinceton", "PRI",
                                           ifelse(method.name=="ASD_frn", "FRN",
                                                  ifelse(method.name=="DAMAGES", "DAM",
                                                         ifelse(method.name=="DAWN", "DAW",
                                                                ifelse(method.name=="RF_Lin", "RFL",
                                                                       ifelse(method.name=="forecASD", "FOR", 
                                                                              ifelse(method.name=="PANDA", "PAN", method.name))))))))%>%
        dplyr::mutate(method.name=ifelse(method.name=="DeRubeis", "DER",
                                         ifelse(method.name=="Sanders", "SAN",
                                                ifelse(method.name=="Satterstrom", "SAT",
                                                       ifelse(method.name=="Iossifov", "IOS", method.name))))) %>%
        dplyr::mutate(method.name = ifelse(method.name == "redo_forec", "Redo",
                                           ifelse(method.name=="noClass_forec", "NoC",
                                                  ifelse(method.name=="noClassPPI_forec","NoCP",
                                                         ifelse(method.name=="noClassPPIBS_forec", "NoCPB",
                                                                ifelse(method.name=="PPIonly_forec", "PPI",
                                                                       ifelse(method.name=="BrainSpanOnly_forec", "BS", method.name))))))) %>%
        dplyr::mutate(method.name= factor(method.name, levels=levels2))%>%
        separate(., AUC_CI, c("p5AUC", "p95AUC"), ", ") %>% 
        separate(., PR20_CI, c("p5PR", "p95PR"), ", ") %>%
        dplyr::mutate(p5AUC = as.numeric(p5AUC), p95AUC=as.numeric(p95AUC))%>%
        dplyr::mutate(p5PR = as.numeric(p5PR),
                      p95PR=as.numeric(p95PR), 
                      realP20R = as.numeric(gsub("\\*$", "", realP20R))) %>%
        #dplyr::mutate(colour = brewer.pal(7, "Dark2"))%>%
        dplyr::select(method.name, realAUC , p5AUC, p95AUC) %>%
        unique%>%
        ggplot(aes(method.name, realAUC))+
        geom_point(aes(colour=method.name), size=1)+
        #facet_wrap(~category,scales="free_x", nrow=1)+
        ylim(0,1)+
        geom_errorbar(aes(ymin=p5AUC, ymax=p95AUC, width=.2),
                      position = position_dodge(0.05), size=1,colour=colours) +
        scale_colour_manual(values = colours) +
        xlab("")+
        ylab(paste(ylabB, "AUROC", sep="\n"))+
        labs(#title = paste0("B) ", title1, "AUROC 95% Confidence Intervals"),
            title = title2, 
            col = "Paper",
            x="Paper")+
        theme_bw()+
        guides(col=guide_legend(nrow=3))+
        theme(legend.title = element_blank())
    #theme(plot.title = element_text(size=34, family="Helvetica"))+
    #theme(axis.text = element_text(size=26, family="Helvetica"),
    #      axis.title = element_text(size=28, family="Helvetica"),
    #      legend.text = element_text(size=24, family="Helvetica"),
    #      legend.title = element_text(size=2, family="Helvetica"),
    #      legend.position="none")
    
    #png = paste0(type, " test.png")
    
    
    
    
    #png(png1, width = 700, height = 1000)
    #gridExtra::grid.arrange(g2,g1,nrow=2)
    #dev.off()
    
    #png(png2, width = 700, height = 1000)
    #gridExtra::grid.arrange(g4,g3,nrow=2)
    #dev.off()
    
    #png(png, width=1400, height=1000)
    g=ggarrange(g2,g4,g1,g3,ncol = 1, nrow = 4, common.legend = TRUE,legend="bottom")
    #dev.off()
    return(g)
    
}

###
#TADA_novel; test expt
###
gba = plotROC_PR(type="GBA_ML",
                 score = c("ASDprinceton_score", "ASD_frn_score", "DAMAGES_score","RF_Lin_score", "forecASD_score", "DAWN_score", "PANDA_score"),
                 genesInPhenotype = rep("TADA_novel", 7),
                 positiveLabels = c("ASDprinceton_posLabs", "ASDfrn_posLabs", "DAMAGES_posLabs", "ASDfrn_posLabs", "forecASD_posLabs", "new.SFARIHC", "PANDA_posLabs"),
                 SFARIset = rep("new.SFARIHC",7),
                 sampleSet= rep("",7))

gen = plotROC_PR(type="Genetics",
                 score = c("DeRubeis_score", "Sanders_score", "iHart_score","Satterstrom_score", "Iossifov_score"),
                 genesInPhenotype = rep("TADA_novel", 5),
                 positiveLabels = c("new.SFARIHC", "new.SFARIHC", "new.SFARIHC", "new.SFARIHC", "new.SFARIHC"),
                 SFARIset = rep("new.SFARIHC",5),
                 sampleSet= rep("",5))

ger = plotROC_PR(type="Generic",
                 score = c("exac_pLI_score", "gnomad_pLI_score", "oe_lof_upper_score"),
                 genesInPhenotype = rep("TADA_novel", 3),
                 positiveLabels = c( "new.SFARIHC", "new.SFARIHC", "new.SFARIHC"),
                 SFARIset = rep("new.SFARIHC",3),
                 sampleSet= rep("",3))
#dev.off()

grDevices::postscript("./Results/Plots/novel_AUROC_PR_CI.ps")
ggarrange(gba, ger, gen, ncol=3)
dev.off()

###
#SFARIHC; cntrl expt
###

gba = plotROC_PR(type="GBA_ML",
                 score = c("ASDprinceton_score", "ASD_frn_score", "DAMAGES_score","RF_Lin_score", "forecASD_score", "DAWN_score", "PANDA_score"),
                 genesInPhenotype = rep("new.SFARIHC", 7),
                 positiveLabels = rep("", 7),
                 SFARIset = rep("",7),
                 sampleSet= rep("",7))

gen = plotROC_PR(type="Genetics",
                 score = c("DeRubeis_score", "Sanders_score", "iHart_score","Satterstrom_score", "Iossifov_score"),
                 genesInPhenotype = rep("new.SFARIHC", 5),
                 positiveLabels = rep("", 5),
                 SFARIset = rep("",5),
                 sampleSet= rep("",5))

ger = plotROC_PR(type="Generic",
                 score = c("exac_pLI_score", "gnomad_pLI_score", "oe_lof_upper_score"),
                 genesInPhenotype = rep("new.SFARIHC", 3),
                 positiveLabels = rep("", 3),
                 SFARIset = rep("",3),
                 sampleSet= rep("",3))

grDevices::postscript("./Results/Plots/SFARIHC_AUROC_PR_CI.ps")
ggarrange(gba, ger, gen, ncol=3) 

dev.off()

###
#forecASD models on TADAnovel; cntrl expt
###
forec_novel = plotROC_PR(type="forec",
                         score = c("redo_forec_score", "noClass_forec_score", "noClassPPI_forec_score","noClassPPIBS_forec_score", "PPIonly_forec_score", "BrainSpanOnly_forec_score"),
                         genesInPhenotype = rep("TADA_novel", 6),
                         positiveLabels = rep("forecASD_posLabs", 6),
                         SFARIset = rep("new.SFARIHC",6),
                         sampleSet= rep("",6))
#dev.off()


###
#forecASD models on SFARIHC; cntrl expt
###
forec_sfari = plotROC_PR(type="forec",
                         score = c("redo_forec_score", "noClass_forec_score", "noClassPPI_forec_score","noClassPPIBS_forec_score", "PPIonly_forec_score", "BrainSpanOnly_forec_score"),
                         genesInPhenotype = rep("new.SFARIHC", 6),
                         positiveLabels = rep("",6),
                         SFARIset = rep("",6),
                         sampleSet= rep("",6))


grDevices::postscript("./Results/Plots/forec_AUROC_PR_CI.ps")
ggarrange(forec_novel, forec_sfari, ncol=2)
dev.off()



