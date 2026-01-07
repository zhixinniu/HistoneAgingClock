############################################
## Load required libraries
############################################
# Description:
# Load all R packages required for downstream ChIP-seq signal processing,
# statistical analysis, visualization, and elastic-net modeling.


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




############################################
## Calculate signal value of each consensus peak
############################################
# Description:
# Quantify ChIP-seq signal for each consensus peak across all samples
# using TCseq. Reads overlapping each peak are counted from BAM files,
# normalized to RPKM, and log-transformed.
#
# Input:
# - final_consensus_peak_set.bed: consensus peak genomic coordinates
# - chip_metadata: metadata table containing sample IP names and Age
# - BAM files: aligned ChIP-seq reads for each sample
#
# Output:
# - tca: TCseq TCA object storing read counts and normalized signal
# - chip_all_signal: matrix of normalized ChIP-seq signal values

# Read consensus peak coordinates
final_consensus_peak_set <- read.table('data/peakcalling/final_consensus_peak_set.bed',header = F)
final_consensus_peak_set$id <- paste0('peak_',1:nrow(final_consensus_peak_set))
colnames(final_consensus_peak_set)[1:3] <- c('chr','start','end')

# Define the directory containing alignment BAM files
path_to_bam <- normalizePath('.data/Alignment/')  # Path to alignment BAM file.

# Construct design table required by TCseq
# - sampleid: ChIP sample identifier
# - timepoint: age encoded as time-course variable
# - group: grouping samples by age
# - BAMfile: corresponding BAM file name
design <- data.frame(sampleid=chip_metadata$IP,
                     timepoint=paste0('year_',chip_metadata$Age),
                     group=chip_metadata%>%
                       group_by(Age)%>%
                       mutate(group_id=paste0('g',cur_group_id()))%>%
                       pull(),
                     BAMfile=paste0(chip_metadata$IP,'.unique_alignment_sorted_rd.bam'))  # Create design information for TCseq.

# Initialize TCseq TCA object linking peaks with sample design
tca <- TCA(design = design,genomicFeature = final_consensus_peak_set)

# Register parallel backend
# Note: number of workers is set to 1 to avoid excessive memory usage
register(MulticoreParam(workers = 1))  # The number of cores used. High RAM usage for large BAM files.

# Count reads overlapping each consensus peak from BAM files
tca <- TCseq::countReads(tca, dir = path_to_bam)  # Count reads for each peak.
tca <- DBanalysis(tca)
tca <- timecourseTable(tca,value = 'expression',filter = F,norm.method = 'rpkm')

# Get counts table for each peak and each sample.
chip_all_signal <- TCseq::counts(tca,normalization='rpkm',log=T)  




############################################
## Remove batch effect
############################################
# Description:
# Correct batch effects (e.g. lab-specific variation) in normalized
# ChIP-seq signal using ComBat, while preserving age-associated effects.
#
# Input:
# - chip_all_signal: normalized ChIP-seq signal matrix
# - metadata.tsv: sample metadata including batch (Lab) and Age
#
# Output:
# - signal_combat_corrected: batch-corrected signal matrix
# - PCA result for evaluating batch effect removal

# Load sample metadata
chip_metadata <- read.table('data/metadata.tsv',sep = '\t',header = T)
batch_metadata <- chip_metadata
rownames(batch_metadata) <- paste0(batch_metadata$IP,'.unique_alignment_sorted.bam')

# Remove non-covariate columns (e.g. sample name, IP)
batch_metadata <- batch_metadata[,-c(1:2),drop=F]

# Apply ComBat to remove batch effects
# - batch: lab information
# - mod: age included as biological covariate
signal_combat_corrected <- ComBat(dat = chip_all_signal,
                                  batch = batch_metadata$Lab,
                                  mod = model.matrix(~batch_metadata$Age,data = meta),
                                  par.prior = T,
                                  prior.plots = F)

# Perform PCA on batch-corrected signal for quality control
pca_res <- PCA(t(signal_combat_corrected),scale.unit = T,graph = F)

fviz_pca_ind(pca_res,
             geom.ind = "point",  
             col.ind = batch_metadata$Lab,  
             palette = "jco",     
             addEllipses = F,  
             legend.title = "Project")




############################################
## Identify age-correlated peaks
############################################
# Description:
# Identify consensus peaks whose ChIP-seq signal changes monotonically
# with age using Spearman correlation.
#
# Input:
# - tca: TCseq TCA object containing normalized time-course signal
#
# Output:
# - chip_cor_peak: table of peaks with correlation coefficient,
#   p-value, and direction of age-associated regulation

# Define function to compute age–signal correlation per peak
cor_peak_ident <- function(tca.obj){

  # Extract time-course signal table from TCseq object
  tcTable <- as.data.frame(tca.obj@tcTable)
  tcTable$id <- rownames(tcTable)
  tcTable <- melt(tcTable,id.vars='id')

  tcTable <- tcTable%>%
    mutate(variable=as.numeric(str_replace(variable,'year_','')))%>%
    mutate(value=as.numeric(value))

  # Compute Spearman correlation for each peak independently
  corTable <- split(tcTable,tcTable$id)%>%
    lapply(.,function(df){
      cor.res <- cor.test(df$value,df$variable,method='spearman',exact = F)
      return(data.frame(peak=df[1,'id'],
                        pval=cor.res$p.value,
                        estimate=unname(cor.res$estimate),
                        regulation=ifelse(cor.res$estimate >= 0.5,'greater',
                                           ifelse(cor.res$estimate <= -0.5,'less','none'))))
    })%>%
    do.call(rbind,.)

  return(corTable)
}

# Calculate age–signal correlations for all consensus peaks
chip_cor_peak <- cor_peak_ident(tca)  




############################################
## Visualization of age-correlated peaks
############################################
# Description:
# Visualize the distribution of Spearman correlation coefficients
# and summarize the proportion of positively, negatively, and
# non-correlated peaks.
#
# Input:
# - chip_cor_peak: age–signal correlation table
#
# Output:
# - Histogram of correlation coefficients with category annotations

# Calculate proportion of peaks in each correlation category
peak_prop <- c(round(sum(chip_cor_peak$alternative=='greater')/nrow(chip_cor_peak),digits = 3)*100,
               round(sum(chip_cor_peak$alternative=='less')/nrow(chip_cor_peak),digits = 3)*100,
               round(sum(chip_cor_peak$alternative=='none')/nrow(chip_cor_peak),digits = 3)*100)
peak_prop <- paste0(peak_prop,'%')

ggplot(chip_cor_peak,aes(x=estimate))+
  geom_histogram(aes(fill=alternative),bins=80,color='white',size=0.05)+
  annotate('label',x=0.75,y=300,label=peak_prop[1],size=5,color='red')+
  annotate('label',x= -0.75,y=300,label=peak_prop[2],size=5,color='blue')+
  annotate('label',x=0,y=300,label=peak_prop[3],size=5,color='gray50')+
  scale_fill_manual(values = c('#F45050','#070A52','gray'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Spearman's r", y='Number of age-correlated peaks')




############################################
## Elastic-net modeling
############################################
# Description:
# Construct an elastic-net regression model to predict chronological age
# based on age-associated ChIP-seq peaks.
#
# Input:
# - chip_all_signal: normalized ChIP-seq signal matrix
# - chip_cor_peak: table of age-correlated peaks
#
# Output:
# - Elastic-net regression model
# - Age predictions and model performance metrics

# Helper function to extract optimal hyperparameters from CV results
get_model_params <- function(fit) {
  alpha <- fit$alpha
  lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
             lambdaSE = lambdaSE[best], eror = error[best])
}  # Output the best hyper-parameters


############################################
## Training/test split, model fitting, and evaluation
############################################
# Description:
# Subset age-associated peaks, split samples into training and test sets,
# perform cross-validation, fit the elastic-net model, and evaluate
# prediction accuracy.
#
# Input:
# - chip_all_signal
# - chip_cor_peak
#
# Output:
# - predictions: observed vs predicted age
# - model.metrics: RMSE, MAE, and Pearson correlation

# Subset signal matrix to age-associated peaks and transpose to samples × peaks
chip_all_signal.sub <- t(chip_all_signal[rownames(chip_all_signal)%in%rownames(subset(chip_cor_peak,alternative!='none')),])

# Generate training dataset
train_index <- data.frame(group=rownames(chip_all_signal.sub),
                          index=1:nrow(chip_all_signal.sub))%>%
  group_by(group)%>%
  sample_frac(size=0.7)%>%
  pull(index)

chip_all_signal.sub.train <- chip_all_signal.sub[train_index,]
rownames(chip_all_signal.sub.train) <- as.numeric(str_replace(rownames(chip_all_signal.sub.train),'year_',''))

# Generate test dataset
chip_all_signal.sub.test <- chip_all_signal.sub[-train_index,]
rownames(chip_all_signal.sub.test) <- as.numeric(str_replace(rownames(chip_all_signal.sub.test),'year_',''))

# Perform cross-validation to select hyperparameters
model.cva <- cva.glmnet(chip_all_signal.sub.train,
                        as.numeric(rownames(chip_all_signal.sub.train)),
                        alpha = seq(0.1,0.9,0.1),
                        nfolds = 10)

# Extract optimal alpha and lambda
hp <- get_model_params(model.cva)

# Fit final elastic-net model
model.glm <- glmnet(chip_all_signal.sub.train,
                    as.numeric(rownames(chip_all_signal.sub.train)),
                    family='gaussian',
                    alpha=hp$alpha)

# Predict age on test dataset
predictions <- as.data.frame(predict(model.glm,
                                     chip_all_signal.sub.test,
                                     type = 'response',
                                     s=hp$lambdaMin))
predictions <- data.frame(real.age=as.numeric(rownames(predictions)),
                          predicted.age=as.numeric(predictions$s1))

# Calculate model performance metrics
model.metrics <- data.frame(RMSE=round(rmse(predictions$real.age,predictions$predicted.age),digits = 3),
                            MAE=round(mae(predictions$real.age,predictions$predicted.age),digits=3),
                            R=round(cor(predictions$real.age,predictions$predicted.age,method = 'pearson'),digits=3))

ggplot(predictions,aes(x=real.age,y=predicted.age))+
  geom_smooth(method = 'lm')+
  geom_point()+
  geom_abline(slope = 1)+
  annotate('label',
           x=min(predictions$real.age),
           y=max(predictions$predicted.age)-5,
           hjust=0,vjust=0,
           label=paste0('RMSE = ',model.metrics$RMSE,
                        '\nMAE = ',model.metrics$MAE,
                        '\nR = ',model.metrics$R))







############################################
## Estimation of sample saturation
############################################
# Description:
# Perform sample-size saturation analysis for histone modification–based
# aging clock models. This analysis evaluates how prediction performance
# (RMSE, MAE, Pearson's r) scales with increasing training sample size.
#
# Input:
# - chip_all_signal.sub.train: original training dataset
# - chip_all_signal.sub.test: independent test dataset
#
# Output:
# - Pearson's r of two independant datasets across sample sizes
# - Observed saturation curve and estimated saturation curve

# Randomly get ID of datasets
id_origin <- nrow(chip_all_signal.sub.train)
id_origin <- sample(id_origin)
id_A <- id_origin[1:floor(length(id_origin)/2)]
id_B <- id_origin[-id_A]

id_origin <- 1:nrow(chip_all_signal.sub.train)
id_doubled <- sample(c(id_origin,id_origin))
id_A_doubled <- id_doubled[1:floor(length(id_doubled)/2)]
id_B_doubled <- id_doubled[-id_A_doubled]


# Calculate Pearson's r of the prediction results of two independant datasets 
saturation_ori <- lapply(seq(10,length(id_A),1), function(x){
  sample_size <- x

  train <- lapply(1:5, function(x){
    train_A <- chip_all_signal.sub.train[sample(id_A,sample_size),]
    train_B <- chip_all_signal.sub.train[sample(id_B,sample_size),]

    return(list(train_A = train_A,
                train_B = train_B))
  })

  print(sample_size)

  res.list <- lapply(train, function(x){
    ## train_A
    cv.fit <- cva.glmnet(x$train_A,as.numeric(rownames(x$train_A)),alpha = seq(0.1,0.9,0.1))
    lambda <- unname(get_model_params(cv.fit)[1,2])
    alpha <- unname(get_model_params(cv.fit)[1,1])
    fit <- glmnet(x$train_A,as.numeric(rownames(x$train_A)),family = 'gaussian',alpha = alpha, nlambda = 100)
    pred_train_A <- predict(fit, chip_all_signal.sub.test, type = 'response',s=lambda)
    pred_train_A <- as.data.frame(pred_train_A)
    pred_train_A <- data.frame('predicted.age'=as.numeric(pred_train_A$s1),
                               'real.age'=as.numeric(blood_k27ac_env$blood_k27ac_test.dat$age))
    ## train_B
    cv.fit <- cva.glmnet(x$train_B,as.numeric(rownames(x$train_B)),alpha = seq(0.1,0.9,0.1))
    lambda <- unname(get_model_params(cv.fit)[1,2])
    alpha <- unname(get_model_params(cv.fit)[1,1])
    fit <- glmnet(x$train_B,as.numeric(rownames(x$train_B)),family = 'gaussian',alpha = alpha, nlambda = 100)
    pred_train_B <- predict(fit, blood_k27ac_env$blood_k27ac_test.rpkm, type = 'response',s=lambda)
    pred_train_B <- as.data.frame(pred_train_B)
    pred_train_B <- data.frame('predicted.age'=as.numeric(pred_train_B$s1),
                               'real.age'=as.numeric(rownames(blood_k27ac_test.rpkm)))
    ## train res
    r_train <- cor(pred_train_A$predicted.age, pred_train_B$predicted.age)

    res.df <- data.frame(group=paste0('prop_',prop),
                         sample_size = c(sample_size),
                         r=c(r_train),
                         type = 'saturation')
    return(res.df)
    })

  res.list <- as.data.frame(rbindlist(res.list))

  return(res.list)

})

saturation_ori.df <- do.call(rbind,saturation_ori)



# Calculate Pearson's r of the prediction results of two doubled independant datasets 
saturation_est <- lapply(seq(10,length(id_A_doubled),20), function(x){
  sample_size <- x

  train <- lapply(1:5, function(x){
    est_A <- chip_all_signal.sub.train[sample(id_A_doubled,sample_size),]
    est_B <- chip_all_signal.sub.train[sample(id_B_doubled,sample_size),]

    return(list(est_A = est_A,
                est_B = est_B))
  })

  print(sample_size)

  res.list <- lapply(train, function(x){
    ## est_A
    cv.fit <- cva.glmnet(x$est_A,as.numeric(rownames(x$est_A)),alpha = seq(0.1,0.9,0.1))
    lambda <- unname(get_model_params(cv.fit)[1,2])
    alpha <- unname(get_model_params(cv.fit)[1,1])
    fit <- glmnet(x$est_A,as.numeric(rownames(x$est_A)),family = 'gaussian',alpha = alpha, nlambda = 100)
    pred_est_A <- predict(fit, chip_all_signal.sub.test, type = 'response',s=lambda)
    pred_est_A <- as.data.frame(pred_est_A)
    pred_est_A <- data.frame('predicted.age'=as.numeric(pred_est_A$s1),
                             'real.age'=as.numeric(blood_k27ac_env$blood_k27ac_test.dat$age))
    ## est_B
    cv.fit <- cva.glmnet(x$est_B,as.numeric(rownames(x$est_B)),alpha = seq(0.1,0.9,0.1))
    lambda <- unname(get_model_params(cv.fit)[1,2])
    alpha <- unname(get_model_params(cv.fit)[1,1])
    fit <- glmnet(x$est_B,as.numeric(rownames(x$est_B)),family = 'gaussian',alpha = alpha, nlambda = 100)
    pred_est_B <- predict(fit, blood_k27ac_env$blood_k27ac_test.rpkm, type = 'response',s=lambda)
    pred_est_B <- as.data.frame(pred_est_B)
    pred_est_B <- data.frame('predicted.age'=as.numeric(pred_est_B$s1),
                             'real.age'=as.numeric(rownames(blood_k27ac_test.rpkm)))
    ## train res
    r_est <- cor(pred_est_A$predicted.age, pred_est_B$predicted.age)

    res.df <- data.frame(group=paste0('prop_',prop),
                         sample_size = c(sample_size),
                         r=c(r_est),
                         type = 'Estimated_saturation')
    return(res.df)
  })

  res.list <- as.data.frame(rbindlist(res.list))

  return(res.list)

})

saturation_est.df <- do.call(rbind,saturation_est)

# Merge final results
saturation_all.df <- rbind(saturation_est.df,saturation_ori.df)

# Plot Pearson's r
ggplot(saturation_all.df, aes(x=sample_size, y=r,color=type, group=type))+
  geom_smooth(span=1,se=F)


