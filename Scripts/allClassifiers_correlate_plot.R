###
#CORRELATE AND HEATMAP all classifiers
###


###
#libraries
###
library(dplyr)
library(ggplot2)
library(zoo)

###
#data and meta
###

allClassifiers_cor = read.csv("Data/allClassifiers/transformedClassifiers/allClassifiers_202012.csv", stringsAsFactors = FALSE) %>%
    dplyr::select(-X) 




###
#correlation heatmap
###
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat, diag=T)]<- NA
    return(cormat)
}


labelColourY = c(ASD_frn_score = "#E41A1C", 
                 DAMAGES_score= "#E41A1C",
                 RF_Lin_score = "#E41A1C",
                 PANDA_score = "#E41A1C",
                 forecASD_score = "#E41A1C", 
                 DAWN_score = "#E41A1C",
                 DeRubeis_score = "#377EB8",
                 Sanders_score = "#377EB8",
                 iHart_score = "#377EB8",
                 Satterstrom_score = "#377EB8",
                 Iossifov_score = "#377EB8",
                 new.SFARIgeneScore = "#4DAF4A",
                 exac_pLI = "#984EA3", 
                 gnomad_pLI = "#984EA3", 
                 oe_lof_upper = "#984EA3", 
                 mis_z = "#984EA3",
                 mfx.Rank = "#FF7F00", 
                 numPubs = "#FF7F00", 
                 physical_node_degree = "#FF7F00")

labelColourX = c(ASDprinceton_score = "#E41A1C",
                 ASD_frn_score = "#E41A1C", 
                 DAMAGES_score= "#E41A1C",
                 RF_Lin_score = "#E41A1C",
                 PANDA_score = "#E41A1C",
                 forecASD_score = "#E41A1C", 
                 DAWN_score = "#E41A1C",
                 DeRubeis_score = "#377EB8",
                 Sanders_score = "#377EB8",
                 iHart_score = "#377EB8",
                 Satterstrom_score = "#377EB8",
                 Iossifov_score = "#377EB8",
                 new.SFARIgeneScore = "#4DAF4A",
                 exac_pLI = "#984EA3", 
                 gnomad_pLI = "#984EA3", 
                 oe_lof_upper = "#984EA3", 
                 mis_z = "#984EA3",
                 mfx.Rank = "#FF7F00", 
                 numPubs = "#FF7F00")


#looking at: 
#all GBA ML
#all hybird
#all GA
#new.SFARIscores
#exac and gnomad and pLI, and oe lof
#N.Lim Mfx scores used
#S.Bhuiyan pubmed numbers used
#biogrid PPI used

#png("Results/allClassifiers/plots/paper_heatmap_all_testFinal_20201016.png", width=825, height=563)
grDevices::postscript("Results/allClassifiers/plots/allScoreHeatmap.ps")
allClassifiers_cor %>%
    dplyr::select(primary.gene.symbol,
                  ASDprinceton_score, ASD_frn_score, DAMAGES_score,RF_Lin_score,PANDA_score,
                  forecASD_score, DAWN_score,
                  DeRubeis_score,Sanders_score,iHart_score,Satterstrom_score,Iossifov_score,
                  new.SFARIgeneScore,
                  exac_pLI_score, gnomad_pLI_score, oe_lof_upper_score, mis_z, mfx.Rank, numPubs, physical_node_degree)%>%
    dplyr::rename_at(vars(contains("_score")), list(~gsub("_score", "", .))) %>%
    dplyr::rename_at(vars(contains("princeton")), list(~gsub("ASDp", "P", .))) %>%
    dplyr::rename_at(vars(contains("ASD_frn")), list(~gsub("ASD_frn", "FRN", .))) %>%
    dplyr::rename_at(vars(contains("new.SFARIgeneScore")), list(~gsub("new.", "", .))) %>%
    dplyr::rename_at(vars(contains("oe_lof_upper")), list(~gsub("oe_lof_upper", "oe_LoF", .))) %>%
    dplyr::rename_at(vars(contains("exac_pLI")), list(~gsub("exac", "ExAC", .))) %>%
    dplyr::rename_at(vars(contains("gnomad_pLI")), list(~gsub("gnomad", "gnomAD", .))) %>%
    dplyr::mutate_at(vars(Princeton:physical_node_degree), ~ ifelse(is.na(.x), na.fill(.x, 0), .x))%>%
   
    #) %>% 
    textshape::column_to_rownames("primary.gene.symbol") %>%
    cor(method="spearman", use='pairwise.complete.obs') %>% 
    get_upper_tri(.) %>%
    reshape2::melt(.,na.rm = TRUE) %>%
    ggplot(., aes(x=Var1, y=Var2, fill=value)) +
    geom_tile(colour="white")+
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint=0, limit=c(-1,1), space="Lab",
                         name="Spearman\nCorrelation")+
    
    #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
    #                             title.position = "top", title.hjust = 0.5))+
    
    
    theme_minimal()+
    theme(axis.text.x=element_text(angle=45, vjust=1,
                                   size=12, hjust=1)) +
    coord_fixed()+
    geom_text(aes(Var1, Var2, label=round(value,2)), colour="black", size=2.5)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y=(element_text(colour=labelColourY, size=12)),
          axis.text.x=(element_text(colour=labelColourX, size=12)),
          legend.justification = c(0, 1),
          legend.position = c(1, 0.9),
          legend.direction = "vertical") +
    annotate("text",  x=17,y=8, label = "Score", fontface =2,size=4,hjust = 0)+
    annotate("text",  x=17,y=7.25, label = c("GBA ML"),size=4,hjust = 0) +
    annotate("text",  x=17,y=6.5, label = c("GA"),size=4,hjust = 0) +
    annotate("text",  x=17,y=5.75, label = c("SFARI"),size=4,hjust = 0) +
    annotate("text",  x=17,y=5, label = c("Constraint"),size=4,hjust = 0) +
    annotate("text",  x=17,y=4.25, label = c("Generic"),size=4,hjust = 0) +
    annotate("text", x=16.5, y=c(7.25,6.5,5.75,5,4.25), label=c("–","–","–","–","–"),fontface=2, size=4, hjust=0) +
    annotate("text",  x=14.8,y=8, label = "Colour", fontface =2,size=4, hjust=0) +
    annotate("rect", xmin = 15, xmax = 16.2, ymin = 7, ymax = 7.5,fill="#E41A1C")+
    annotate("rect", xmin = 15, xmax = 16.2, ymin = 6.25, ymax = 6.75,fill="#377EB8") +
    annotate("rect", xmin = 15, xmax = 16.2, ymin = 5.5, ymax = 6,fill="#4DAF4A") +
    annotate("rect", xmin = 15, xmax = 16.2, ymin = 4.75, ymax = 5.25,fill="#984EA3") +
    annotate("rect", xmin = 15, xmax = 16.2, ymin = 4, ymax = 4.5,fill="#FF7F00") 
    

dev.off()



###
#overlap
###
allClassifiers_overlap = allClassifiers_cor %>%
    dplyr::select(primary.gene.symbol, ASDprinceton_top, ASD_frn_top, DAMAGES_top, RF_lin_top, forecASD_top, DAWN_top, 
                  DeRubeis_top, Sanders_top, iHart_top, Satterstrom_top, Iossifov_top,
                  exac_pLI_top, gnomad_pLI_top, oe_lof_upper_top, new.SFARIHC, TADA_novel) %>%
    dplyr::rename(., SFARIHC = new.SFARIHC) %>%
    dplyr::rename(.,  NOVELHC = TADA_novel) 


overlap = crossprod(as.matrix(allClassifiers_overlap[2:17]))            # calculate the overlap               
overlap[lower.tri(overlap)] <- NA

write.csv(overlap, "Results/allClassifiers/topOverlap.csv")


