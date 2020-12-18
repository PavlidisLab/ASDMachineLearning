#### load packages ####
# # # # # # # # # # # # #
#### Ensembel Model #####
# # # # # # # # # # # # #

library(randomForest)
library(dplyr)

unzip("./Data/forecASD_data.zip")

load("./Data/forecASD_data/01_training_labels.Rdata")
load("./Data/forecASD_data/02_STRING_rf.Rdata")
load("./Data/forecASD_data/02_brainspan_RF.Rdata")
load("./Data/forecASD_data/02_network_rf.Rdata")

# # # # # # # # # # # # 
#### Load in data #####
# # # # # # # # # # # # 

meta = read.csv("./Data/forecASD_data/composite_table.csv", stringsAsFactors = F, row.names = 1)

string.prd = string.prd[rownames(meta), ]
bs.prd = bs.prd[rownames(meta), ]

## combine other predictors with network scores
meta = cbind(
    data.frame(
        STRING_score = string.prd[rownames(meta) , "TRUE"],
        BrainSpan_score = bs.prd[rownames(meta), "TRUE"]
    ),
    meta[rownames(meta), ]
)

names(meta)[14:17]=c("tada_asc+ssc+del","tada_asc+ssc", "tada_asc", "tada_ssc")
names(meta)[12:13]=c("D_ens","D")
# # # # # # # # # # # # 
#### train forest #####
# # # # # # # # # # # # 

meta.train = na.roughfix(
    meta[meta$ensembl_string %in% c(pos,neg), -(3:9)]   ## remove gene identifiers, etc.
) 
y = as.factor(rownames(meta.train) %in% pos)

set.seed(43775)
rf1 = randomForest(
    y = y,
    x = meta.train, 
    importance = T,
    do.trace = 10,
    strata = y,
    sampsize = c(76,76)
)

meta.test = na.roughfix(
    meta[!(meta$ensembl_string %in% c(pos,neg)), -(3:9)]   ## remove gene identifiers, etc.
) 

meta.prd <- predict(rf1, 
                    meta.test, 
                    type = "prob")

meta.score <- rbind(rf1$votes, meta.prd)

final.data <- cbind(
    data.frame(
        forecASD = meta.score[rownames(meta),"TRUE"],
        'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
    ), 
    meta
)

final.data = final.data %>% dplyr::select(forecASD, entrez, symbol) %>% dplyr::rename_at(vars(forecASD), list(~gsub("forecASD", "fullForec_score",.)))

# # # # # # # # # # # # 
#### train forest #####: no classifiers
# # # # # # # # # # # # 

meta.train2 = meta.train %>% dplyr::select(-c(krishnan_post,D,D_ens,rASD,min_lFDR,Netscore))

y2 = as.factor(rownames(meta.train2) %in% pos)

set.seed(43775)
rf2 = randomForest(
    y = y2,
    x = meta.train2, 
    importance = T,
    do.trace = 10,
    strata = y2,
    sampsize = c(76,76)
)

meta.test2 = meta.test %>% dplyr::select(-c(krishnan_post,D,D_ens,rASD,min_lFDR,Netscore))


meta.prd <- predict(rf2, 
					meta.test2, 
					type = "prob")

meta.score <- rbind(rf2$votes, meta.prd)

final.data2 <- cbind(
	data.frame(
		forecASD = meta.score[rownames(meta),"TRUE"],
		'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
	), 
	meta
)


final.data2= final.data2 %>% dplyr::select(forecASD, entrez, symbol) %>% dplyr::rename_at(vars(forecASD), list(~gsub("forecASD", "noClass_Forec_score",.)))


# # # # # # # # # # # # 
#### train forest #####: no classifiers, PPI
# # # # # # # # # # # # 

meta.train6 = meta.train %>% dplyr::select(-c(krishnan_post,D,D_ens,rASD,min_lFDR,Netscore, STRING_score))

y6 = as.factor(rownames(meta.train6) %in% pos)

set.seed(43775)
rf6 = randomForest(
	y = y6,
	x = meta.train6, 
	importance = T,
	do.trace = 10,
	strata = y6,
	sampsize = c(76,76)
)

meta.test6 = meta.test %>% dplyr::select(-c( krishnan_post,D,D_ens,rASD,min_lFDR,Netscore, STRING_score))


meta.prd <- predict(rf6, 
					meta.test6, 
					type = "prob")

meta.score <- rbind(rf6$votes, meta.prd)

final.data6 <- cbind(
	data.frame(
		forecASD = meta.score[rownames(meta),"TRUE"],
		'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
	), 
	meta
)


final.data6= final.data6 %>% dplyr::select(forecASD, entrez, symbol) %>% dplyr::rename_at(vars(forecASD), list(~gsub("forecASD", "noPPIClass_Forec_score",.)))


# # # # # # # # # # # # 
#### train forest #####: PPIonly
# # # # # # # # # # # # 

meta.train7 = meta.train %>% dplyr::select(STRING_score)

y7 = as.factor(rownames(meta.train7) %in% pos)

set.seed(43775)
rf7 = randomForest(
    y = y7,
    x = meta.train7, 
    importance = T,
    do.trace = 10,
    strata = y7,
    sampsize = c(76,76)
)

meta.test7 = meta.test %>% dplyr::select(STRING_score)


meta.prd <- predict(rf7, 
                    meta.test7, 
                    type = "prob")

meta.score <- rbind(rf7$votes, meta.prd)

final.data7 <- cbind(
    data.frame(
        forecASD = meta.score[rownames(meta),"TRUE"],
        'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
    ), 
    meta
)


final.data7= final.data7 %>% dplyr::select(forecASD, entrez, symbol) %>% dplyr::rename_at(vars(forecASD), list(~gsub("forecASD", "PPI_Forec_score",.)))

# # # # # # # # # # # # 
#### train forest #####: only genetics 
# # # # # # # # # # # # 

meta.train8 = meta.train %>% dplyr::select(-c(krishnan_post,D,D_ens,rASD,min_lFDR,Netscore, STRING_score, BrainSpan_score))

y8 = as.factor(rownames(meta.train6) %in% pos)

set.seed(43775)
rf8 = randomForest(
    y = y8,
    x = meta.train8, 
    importance = T,
    do.trace = 10,
    strata = y8,
    sampsize = c(76,76)
)



meta.test8 = meta.test %>% dplyr::select(-c( krishnan_post,D,D_ens,rASD,min_lFDR,Netscore, STRING_score, BrainSpan_score))


meta.prd <- predict(rf8, 
                    meta.test8, 
                    type = "prob")

meta.score <- rbind(rf8$votes, meta.prd)

final.data8 <- cbind(
    data.frame(
        forecASD = meta.score[rownames(meta),"TRUE"],
        'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
    ), 
    meta
)


final.data8= final.data8 %>% dplyr::select(forecASD, entrez, symbol) %>% dplyr::rename_at(vars(forecASD), list(~gsub("forecASD", "onlyGenetics_Forec_score",.)))


# # # # # # # # # # # # 
#### train forest #####: brainspan
# # # # # # # # # # # # 
meta.train9 = meta.train %>% dplyr::select( BrainSpan_score)

y9 = as.factor(rownames(meta.train6) %in% pos)

set.seed(43775)
rf9 = randomForest(
    y = y9,
    x = meta.train9, 
    importance = T,
    do.trace = 10,
    strata = y8,
    sampsize = c(76,76)
)






meta.test9 = meta.test %>% dplyr::select( BrainSpan_score)


meta.prd <- predict(rf9, 
                    meta.test9, 
                    type = "prob")

meta.score <- rbind(rf9$votes, meta.prd)

final.data9<- cbind(
    data.frame(
        forecASD = meta.score[rownames(meta),"TRUE"],
        'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
    ), 
    meta
)


final.data9= final.data9 %>% dplyr::select(forecASD, entrez, symbol) %>% dplyr::rename_at(vars(forecASD), list(~gsub("forecASD", "onlyBrainSpan_Forec_score",.)))



finalDat = left_join(final.data, final.data2) %>% 
    left_join(., final.data6)%>% 
    left_join(., final.data7) %>%
    left_join(., final.data8) %>%
    left_join(., final.data9)


#write.csv(finalDat, "results/rf_4.6.14_final/adaptations/finalForecASD_permuataions_rf_4.6.14_redo_20200401.csv")

# # # # # # # # # # # # 
#### var imp plots #####
# # # # # # # # # # # # 



grDevices::postscript("./Results/Plots/forecASD_varIMP_202012.ps")
par(mfrow=c(2,4))
varImpPlot(rf1, main="", type=1)
varImpPlot(rf1, main="", type=2)
mtext("forecASD/Redo", font=2, cex=1, side = 3, line =-3, adj=0.25, outer = TRUE)
varImpPlot(rf2, main="", type=1)
varImpPlot(rf2, main="", type=2)
mtext("noClass", font=2, cex=1, side = 3, line =-3, adj=0.80, outer = TRUE)
varImpPlot(rf6,  main="", type=1)
varImpPlot(rf6, main="", type=1)
mtext("noClassPPI", font=2, cex=1, side = 3, line =-32, adj=0.25, outer = TRUE)
varImpPlot(rf8, main="", type=1)
varImpPlot(rf8, main="", type=2)
mtext("noClassPPIBS", font=2, cex=1, side = 3, line =-32, adj=0.85, outer = TRUE)
dev.off()



# # # # # # # # # # # # 
#### clean data #####
# # # # # # # # # # # # 

###
#allclassifiers + gene alias info
###

unzip("./Data/allClassifiers_202012.zip")
allClassifiers = read.csv("./Data/allClassifiers_202012.csv", stringsAsFactors = FALSE) %>%
	dplyr::select(-X) 


unzip("./Data/geneAliasInfo_all_2012.zip")
geneAliasInfo = read.csv("./Data/geneAliasInfo_all_202012.csv")

###
#forecASD
###
newForecASDEDIT = finalDat %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(primary.gene.symbol, primary.gene.id, primary.gene.biotype,
				  fullForec_score, noClass_Forec_score, 
				  noPPIClass_Forec_score, PPI_Forec_score, onlyGenetics_Forec_score,onlyBrainSpan_Forec_score)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(fullForec_score = mean(fullForec_score),
					 noClass_Forec_score = mean(noClass_Forec_score),
					 noPPIClass_Forec_score = mean(noPPIClass_Forec_score),
					 PPI_Forec_score = mean(PPI_Forec_score),
					 onlyGenetics_Forec_score = mean(onlyGenetics_Forec_score),
					 onlyBrainSpan_Forec_score = mean(onlyBrainSpan_Forec_score)) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(primary.gene.id = as.character(primary.gene.id))#%>%dim #17879
#dplyr::filter(!(grepl(",", primary.gene.symbol)))%>% dim  #17879

newForecASDEDIT$primary.gene.symbol[duplicated(newForecASDEDIT$primary.gene.symbol)]
newForecASDEDIT$primary.gene.id[duplicated(newForecASDEDIT$primary.gene.id)]

dim(newForecASDEDIT)#17879

redoTop = finalDat[,c("symbol", "entrez", "fullForec_score")] %>% top_n(1787) %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::select(fullForec_score,gene.symbol, ncbi.id, mapping,primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(fullForec_score, primary.gene.symbol, primary.gene.id, primary.gene.biotype)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(fullForec_score = mean(fullForec_score)) %>%
	dplyr::ungroup()


dim(redoTop)#1803

noClassTop = finalDat[,c("symbol", "entrez", "noClass_Forec_score")] %>% top_n(1787) %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::select(noClass_Forec_score,gene.symbol, ncbi.id, mapping,primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(noClass_Forec_score, primary.gene.symbol, primary.gene.id, primary.gene.biotype)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(noClass_Forec_score = mean(noClass_Forec_score)) %>%
	dplyr::ungroup()


dim(noClassTop)#1787

noClassPPITop = finalDat[,c("symbol", "entrez", "noPPIClass_Forec_score")] %>% top_n(1787) %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::select(noPPIClass_Forec_score,gene.symbol, ncbi.id, mapping,primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(noPPIClass_Forec_score, primary.gene.symbol, primary.gene.id, primary.gene.biotype)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(noPPIClass_Forec_score = mean(noPPIClass_Forec_score)) %>%
	dplyr::ungroup()


dim(noClassPPITop)#1785

noClassPPIBSTop = finalDat[,c("symbol", "entrez", "onlyGenetics_Forec_score")] %>% top_n(1787) %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::select(onlyGenetics_Forec_score,gene.symbol, ncbi.id, mapping,primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(onlyGenetics_Forec_score, primary.gene.symbol, primary.gene.id, primary.gene.biotype)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(onlyGenetics_Forec_score = mean(onlyGenetics_Forec_score)) %>%
	dplyr::ungroup()


dim(noClassPPIBSTop)#1785


onlyPPITop = finalDat[,c("symbol", "entrez", "PPI_Forec_score")] %>% top_n(1787) %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::select(PPI_Forec_score,gene.symbol, ncbi.id, mapping,primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(PPI_Forec_score, primary.gene.symbol, primary.gene.id, primary.gene.biotype)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(PPI_Forec_score = mean(PPI_Forec_score)) %>%
	dplyr::ungroup()


dim(onlyPPITop)#1784

onlyBSTop = finalDat[,c("symbol", "entrez", "onlyBrainSpan_Forec_score")] %>% top_n(1787) %>%
	dplyr::mutate(gene.symbol=toupper(symbol), ncbi.id = as.integer(entrez)) %>% 
	dplyr::mutate(gene.symbol = ifelse(is.na(gene.symbol), "", gene.symbol))%>%
	left_join(., geneAliasInfo, by=c("gene.symbol", "ncbi.id")) %>% 
	dplyr::select(onlyBrainSpan_Forec_score,gene.symbol, ncbi.id, mapping,primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::filter(mapping=="single") %>%
	dplyr::filter(primary.gene.biotype%in%c("protein-coding")) %>% #dim #17901
	dplyr::select(onlyBrainSpan_Forec_score, primary.gene.symbol, primary.gene.id, primary.gene.biotype)%>%
	unique %>%
	dplyr::group_by(primary.gene.symbol, primary.gene.id, primary.gene.biotype) %>%
	dplyr::summarise(onlyBrainSpan_Forec_score = mean(onlyBrainSpan_Forec_score)) %>%
	dplyr::ungroup()


dim(onlyPPITop)#1784



allClassifiers = allClassifiers %>%
	dplyr::mutate(primary.gene.id = as.character(primary.gene.id)) %>%
	left_join(., newForecASDEDIT) %>%
	dplyr::mutate(redo_forec_score = ifelse(is.na(fullForec_score), 0, fullForec_score),
				  noClass_forec_score = ifelse(is.na(noClass_Forec_score), 0, noClass_Forec_score),
				  noClassPPI_forec_score = ifelse(is.na(noPPIClass_Forec_score), 0, noPPIClass_Forec_score),
				  noClassPPIBS_forec_score = ifelse(is.na(onlyGenetics_Forec_score),0,onlyGenetics_Forec_score),
				  PPIonly_forec_score = ifelse(is.na(PPI_Forec_score), 0, PPI_Forec_score),
				  BrainSpanOnly_forec_score = ifelse(is.na(onlyBrainSpan_Forec_score), 0, onlyBrainSpan_Forec_score)) %>%
	dplyr::mutate(redoTop = ifelse(primary.gene.symbol %in% redoTop$primary.gene.symbol, 1, 0),
				  noClassTop  = ifelse(primary.gene.symbol %in% noClassTop$primary.gene.symbol, 1, 0),
				  noClassPPITop  = ifelse(primary.gene.symbol %in% noClassPPITop$primary.gene.symbol, 1, 0),
				  noClassPPIBSTop  = ifelse(primary.gene.symbol %in% noClassPPIBSTop$primary.gene.symbol, 1, 0),
				  onlyPPITop  = ifelse(primary.gene.symbol %in% onlyPPITop$primary.gene.symbol, 1, 0),
				  onlyBSTop  = ifelse(primary.gene.symbol %in% onlyBSTop$primary.gene.symbol, 1, 0))

write.csv(allClassifiers, "./Data/allClassifiers_forecASD_202012.csv")






