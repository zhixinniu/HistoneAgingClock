library(jsonlite)
library(dplyr)
library(stringr)
library(data.table)
library(TCseq)
library(BiocParallel)
library(reshape2)
library(ggplot2)
library(glmnet)
library(glmnetUtils)
library(Metrics)
library(sva)
library(factoextra)
library(FactoMineR)




##### Calculate signal value of each consensus peak
# Input:
# consensus peakset
# metadata for TCseq
# Ouput:
# TCseq TCA object with calculated and normalized counts table

final_consensus_peak_set <- read.table('data/peakcalling/final_consensus_peak_set.bed',header = F)
final_consensus_peak_set$id <- paste0('peak_',1:nrow(final_consensus_peak_set))
colnames(final_consensus_peak_set)[1:3] <- c('chr','start','end')
path_to_bam <- normalizePath('.data/Alignment/')  # Path to alignment BAM file.
design <- data.frame(sampleid=chip_metadata$IP,
                     timepoint=paste0('year_',chip_metadata$Age),
                     group=chip_metadata%>%
                       group_by(Age)%>%
                       mutate(group_id=paste0('g',cur_group_id()))%>%
                       pull(),
                     BAMfile=paste0(chip_metadata$IP,'.unique_alignment_sorted_rd.bam'))  # Create design information for TCseq.
tca <- TCA(design = design,genomicFeature = final_consensus_peak_set)
register(MulticoreParam(workers = 1))  # The number of cores used. High RAM usage for large BAM files.
tca <- TCseq::countReads(tca, dir = path_to_bam)  # Count reads for each peak.
tca <- DBanalysis(tca)
tca <- timecourseTable(tca,value = 'expression',filter = F,norm.method = 'rpkm')
chip_all_signal <- TCseq::counts(tca,normalization='rpkm',log=T)  # Get counts table for each peak and each sample.





##### Remove batch effect
# Input:
# normalized signal value
# Output:
# Batch effect-corrected signal value

chip_metadata <- read.table('data/metadata.tsv',sep = '\t',header = T)
batch_metadata <- chip_metadata
rownames(batch_metadata) <- paste0(batch_metadata$IP,'.unique_alignment_sorted.bam')
batch_metadata <- batch_metadata[,-c(1:2),drop=F]


signal_combat_corrected <- ComBat(dat = chip_all_signal,
                                  batch = batch_metadata$Lab,
                                  mod = model.matrix(~batch_metadata$Age,data = meta),
                                  par.prior = T,
                                  prior.plots = F)

pca_res <- PCA(t(signal_combat_corrected),scale.unit = T,graph = F)
fviz_pca_ind(pca_res,
             geom.ind = "point",  
             col.ind = batch_metadata$Lab,  
             palette = "jco",     
             addEllipses = F,  
             legend.title = "Project")




##### Identify age-correlated peaks
# Input:
# signal value table
# Output
# table of potentially age-associated peaks
cor_peak_ident <- function(tca.obj){
  tcTable <- as.data.frame(tca.obj@tcTable)
  tcTable$id <- rownames(tcTable)
  tcTable <- melt(tcTable,id.vars='id')
  tcTable <- tcTable%>%
    mutate(variable=as.numeric(str_replace(variable,'year_','')))%>%
    mutate(value=as.numeric(value))
  corTable <- split(tcTable,tcTable$id)%>%
    lapply(.,function(df){
      cor.res <- cor.test(df$value,df$variable,method='spearman',exact = F)
      return(data.frame(peak=df[1,'id'],
                        pval=cor.res$p.value,
                        estimate=unname(cor.res$estimate),
                        regulation=ifelse(cor.res$estimate >= 0.5 'greater',
                                           ifelse(cor.res$estimate <= -0.5,'less','none'))))
    })%>%
    do.call(rbind,.)
  return(corTable)
}
chip_cor_peak <- cor_peak_ident(tca)  # Calculate Spearman's correlation during aging.

peak_prop <- c(round(sum(chip_cor_peak$alternative=='greater')/nrow(chip_cor_peak),digits = 3)*100,
               round(sum(chip_cor_peak$alternative=='less')/nrow(chip_cor_peak),digits = 3)*100,
               round(sum(chip_cor_peak$alternative=='none')/nrow(chip_cor_peak),digits = 3)*100)  # The proportion of positively- and negatively correlated peaks.
peak_prop <- paste0(peak_prop,'%')
ggplot(chip_cor_peak,aes(x=estimate))+
  geom_histogram(aes(fill=alternative),bins=80,color='white',size=0.05)+
  annotate('label',x=0.75,y=300,label=peak_prop[1],size=5,color='red')+
  annotate('label',x= -0.75,y=300,label=peak_prop[2],size=5,color='blue')+
  annotate('label',x=0,y=300,label=peak_prop[3],size=5,color='gray50')+
  scale_fill_manual(values = c('#F45050','#070A52','gray'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Spearman's r", y='Number of age-correlated peaks') # Visualize the distribution of age-correlated peaks.




##### Modeling
# Input:
# corrected signal table
# Output:
# elastic-net regression model

# Define training and test datasets, e.g. 70% for training and 30% for test.
get_model_params <- function(fit) {
  alpha <- fit$alpha
  lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
             lambdaSE = lambdaSE[best], eror = error[best])
}  # Output the best hyper-parameters

chip_all_signal.sub <- t(chip_all_signal[rownames(chip_all_signal)%in%rownames(subset(chip_cor_peak,alternative!='none')),])
train_index <- data.frame(group=rownames(chip_all_signal.sub),
                          index=1:nrow(chip_all_signal.sub))%>%
  group_by(group)%>%
  sample_frac(size=0.7)%>%
  pull(index)
chip_all_signal.sub.train <- chip_all_signal.sub[train_index,]  # Generate training dataset
rownames(chip_all_signal.sub.train) <- as.numeric(str_replace(rownames(chip_all_signal.sub.train),'year_',''))
chip_all_signal.sub.test <- chip_all_signal.sub[-train_index,]  # Generate test dataset
rownames(chip_all_signal.sub.test) <- as.numeric(str_replace(rownames(chip_all_signal.sub.test),'year_',''))
model.cva <- cva.glmnet(chip_all_signal.sub.train,as.numeric(rownames(chip_all_signal.sub.train)),alpha = seq(0.1,0.9,0.1),nfolds = 10)  # Cross-validation for model evaluation. Nfolds depend on the number of samples.
hp <- get_model_params(model.cva)  # Get the best hyper-parameters
model.glm <- glmnet(chip_all_signal.sub.train,as.numeric(rownames(chip_all_signal.sub.train)),family='gaussian',alpha=hp$alpha)  # Establish the elastic net regression model
predictions <- as.data.frame(predict(model.glm,chip_all_signal.sub.test,type = 'response',s=hp$lambdaMin))
predictions <- data.frame(real.age=as.numeric(rownames(predictions)),
                          predicted.age=as.numeric(predictions$s1))  # Final predictions.
model.metrics <- data.frame(RMSE=round(rmse(predictions$real.age,predictions$predicted.age),digits = 3),
                            MAE=round(mae(predictions$real.age,predictions$predicted.age),digits=3),
                            R=round(cor(predictions$real.age,predictions$predicted.age,method = 'pearson'),digits=3))  # Calculate RMSE, MAE, and Pearson's R
ggplot(predictions,aes(x=real.age,y=predicted.age))+
  geom_smooth(method = 'lm')+
  geom_point()+
  geom_abline(slope = 1)+
  annotate('label',x=min(predictions$real.age),y=max(predictions$predicted.age)-5,hjust=0,vjust=0,
           label=paste0('RMSE = ',model.metrics$RMSE,'\nMAE = ',model.metrics$MAE,'\nR = ',model.metrics$R)) # Visualization of predictions






##### Estimation of sample saturation
# This script performs a sample-size saturation analysis for histone modification–based aging clock models.
# Inspired by the MEDIPS saturation framework, it evaluates how model prediction accuracy (RMSE, MAE, Pearson's r)
# scales with increasing sample size and assesses performance stabiligy using sampling and artificial duplication

# Input: 
# training dataset
# Output:
# RMSE (or MAE, Pearson's r) vs. sample size curve
# Estimated saturationg curve


saturation_sample_size <- lapply(seq(0.1,1,0.01), function(i){       # evaluate performance metrics for sub-samp
  prop <- i # proportion of training datasets
  nit <- # iteration number of sampling
  
  train <- lapply(1:nit, function(x){
    idx <- sample(seq_len(nrow(chip_all_signal.sub.train)), size=floor(prop*nrow(chip_all_signal.sub.train)))
    chip_all_signal.sub.train[idx,]
  })
  
  print(paste0('loop of ',prop))
  
  res.list <- lapply(train, function(x){
    cv.fit <- cva.glmnet(x,as.numeric(rownames(x)),alpha = seq(0.1,0.9,0.1))
    lambda <- unname(hp[1,2])
    alpha <- unname(hp[1,1])
    
    fit <- glmnet(x,as.numeric(rownames(x)),family = 'gaussian',alpha = alpha, nlambda = 100)
    pred <- predict(fit, chip_all_signal.sub.test, type = 'response',s=lambda)
    pred <- as.data.frame(pred)
    pred <- data.frame('predicted.age'=as.numeric(pred$s1),
                       'real.age'=as.numeric(rownames(pred)))
    mae <- mae(pred$real.age, pred$predicted.age)
    rmse <- rmse(pred$real.age, pred$predicted.age)
    r <- cor(pred$real.age, pred$predicted.age,method = 'pearson')
    
    res.df <- data.frame(group=paste0('prop_',prop),
                         prop=prop,
                         r=r,
                         rmse=rmse,
                         mae=mae)=mae)
    return(res.df)
    
  })
  
  res.list <- as.data.frame(rbindlist(res.list))return(res.list)
})

saturation_sample_size.df <- do.call(rbind, saturation_sample_size)





