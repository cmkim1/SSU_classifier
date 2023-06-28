---
  title: "term_project_230612"
output: html_document
author: csbl_cmkim
date: "2023-06-12"
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

##### 0-1. loading library
```{r library}

#install.packages("https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.19.tar.gz", repos = NULL, type="source")
#BiocManager::install("DNAStringSet")
library(dplyr)
library(ggplot2)
library(ISLR)
library(MASS)
library(glmnet)
library(randomForest)
library(gbm)
library(rpart)
library(boot)
library(data.table)
library(ROCR)
library(gridExtra)
library(tidyr)
library(stringr)
library(Biostrings)
library(psych)
library(tidyverse)
library(tibble)
```


##### 0-2. defining function
```{r function}

#binomial_deviance function
binomial_deviance <- function(y_obs, yhat){
  epsilon = 0.0001
  yhat = ifelse(yhat < epsilon, epsilon, yhat)
  yhat = ifelse(yhat > 1-epsilon, 1-epsilon, yhat)
  a = ifelse(y_obs==0, 0, y_obs * log(y_obs/yhat))
  b = ifelse(y_obs==1, 0, (1-y_obs) * log((1-y_obs)/(1-yhat)))
  return(2*sum(a + b))
}

#predict.cv.glmnet function
predict.cv.glmnet=function(object,newx,s=c("lambda.1se","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else
    if(is.character(s)){
      s=match.arg(s)
      lambda=object[[s]]
      names(lambda)=s
    }
  else stop("Invalid form for s")
  predict(object$glmnet.fit,newx,s=lambda,...)
}


#panel.cor function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
} 
```

##### 1-1. loading taxa
```{r taxa}
seq_itgdb_taxa <- read.csv("/Users/chungminkim/Downloads/drive-download-20230604T125759Z-001/seq_itgdb_taxa.txt",
                           sep='\t', header=F)
taxa <- separate(seq_itgdb_taxa, V2, into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ";")
colnames(taxa)[c(1)] <- c('id')
glimpse(taxa)
```

##### 1-2. loading sequences
```{r sequences}
seq_itgdb_seq <- readDNAStringSet("/Users/chungminkim/Downloads/drive-download-20230604T125759Z-001/seq_itgdb_seq.fasta")
id_sequence = data.frame(id=names(seq_itgdb_seq),sequences=as.character(seq_itgdb_seq))
glimpse(id_sequence)

# 486640 16S sequence in total

```

##### 1-3. processing data & counting frequency
```{r processing}
phylum_seq <- cbind(taxa[,3], id_sequence[,2])
colnames(phylum_seq) <- c('phylum', 'sequence')
phylum_seq <- as_tibble(phylum_seq)
glimpse(phylum_seq)


phylum_seq <- phylum_seq %>% filter(!is.na(phylum))
phylum_seq %>% summarize(n_phylum = n_distinct(phylum)) # number of phylum: 317
taxa %>% group_by(phylum) %>% count(sort=T) # frequency table

phyla = taxa %>% group_by(phylum) %>% count(sort=T)
ggplot(phyla[1:10,], aes(reorder(phylum, -n),n)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# proteobacteria, firmicutes, bacteroidota, actinobacteriota, ....


phylum_table <- phylum_seq %>% group_by(phylum) %>%
  summarize(Freq=n_distinct(sequence))
phylum_table <- phylum_table %>% filter(row_number() <= n()-1) # remove NA


# number of phyla which have more than 100, 1000, 10000 freq
phylum_table %>%
  summarize(more_than_100 = n_distinct(phylum_table %>% filter(Freq > 100)),
            more_than_1000 = n_distinct(phylum_table %>% filter(Freq > 1000)),
            more_than_10000 = n_distinct(phylum_table %>% filter(Freq > 10000))
  ) #120, 38, 5
```


##### 2-1. example1: firmicutes
```{r example1_data}

# randomly selecting 20,000 sequence
data <- phylum_seq %>% sample_n(20000)

#number of firmicutes sequence and others
data %>%
  summarize(n_firmicutes = n_distinct(data %>% filter(phylum == "firmicutes")),
            n_others = n_distinct(data %>% filter(phylum != "firmicutes")))

# counting nucleotide mono, di, trimer
r16Sstr <- DNAStringSet(data$sequence)
k1fq = alphabetFrequency(r16Sstr)
k1fq = k1fq[,1:4]
k2fq = dinucleotideFrequency(r16Sstr)
k3fq = trinucleotideFrequency(r16Sstr)
#k4fq = oligonucleotideFrequency(r16Sstr, width=4)
#k5fq = oligonucleotideFrequency(r16Sstr, width=5)
data <- cbind(data[,1], k1fq, k2fq, k3fq)
colnames(data)[1]<-'phylum'
as_tibble(data)
data$phylum <- factor(ifelse(data$phylum == 'firmicutes', 0, 1))
head(data)
```



##### 2-2. example2 running
```{r running1}
data2 <- data
#data2 <- cbind(data[,1], data_new)
colnames(data2)[1]<-'phylum'

#data.frame(data=dim(data)[2], data2=dim(data2)[2])


#TRAINING
set.seed(2306)
n <- nrow(data2)
idx <- 1:n
training_idx <- sample(idx, n * .60)
idx <- setdiff(idx, training_idx)
validate_idx <- sample(idx, n * .20)
test_idx <- setdiff(idx, validate_idx)
training <- data2[training_idx,]
validation <- data2[validate_idx,]
test <- data2[test_idx,]


#LOGISTIC REGRESSION
data2_lm_full <- glm(phylum ~ ., data=training, family=binomial)
summary.glm(data2_lm_full)
contribution <- as.data.frame(data2_lm_full$coefficients)
contribution %>% rownames_to_column() %>%
  top_n(6, contribution$data2_lm_full.coefficients) %>% pull(rowname)



predict(data2_lm_full, newdata = data2[1:5,], type='response')

#MODEL EVALUATION
y_obs <- as.numeric(as.character(validation$phylum))
yhat_lm <- predict(data2_lm_full, newdata = validation, type='response')
pred_lm <- prediction(as.numeric(yhat_lm), as.numeric(y_obs))
#performance(pred_lm, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_lm) 
#plot(performance(pred_lm, measure="tpr", x.measure='fpr'))


#LASSO MODEL
xx <- model.matrix(phylum ~ .-1, data2)
x <- xx[training_idx, ]
y <- as.numeric(as.character(training$phylum))
#glimpse(x)

library(tictoc)
tic("LASSO")
data2_cvfit <- cv.glmnet(x, y, family = "binomial")
#plot(data2_cvfit)
LASSO_time <- toc()

#coef(data2_cvfit, s = c("lambda.1se"))
#coef(data2_cvfit, s = c("lambda.min")) #예측에는 lambda.min 사용

#MODEL EVALUATION
#predict.cv.glmnet(data2_cvfit, s='lambda.min', newx = x[1:5,], type='response')

yhat_glmnet <- predict(data2_cvfit, s="lambda.min", newx=xx[validate_idx,],
                       type='response')
yhat_glmnet <- yhat_glmnet[,1]
pred_glmnet <- prediction(yhat_glmnet, y_obs)
#performance(pred_glmnet, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_glmnet)

#TREE MODEL
data2_tr <- rpart(phylum ~ ., data = training)
#data2_tr

#printcp(data2_tr)
#summary(data2_tr)

opar <- par(mfrow = c(1,1), xpd = NA)
plot(data2_tr)
text(data2_tr, use.n = T)
par(opar)

#EVALUATION
yhat_tr <- predict(data2_tr, validation)
yhat_tr <- yhat_tr[,"1"]
pred_tr <- prediction(yhat_tr, y_obs)
#performance(pred_tr, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_tr)

#RANDOM FOREST
set.seed(2305)
tic("RF")
data2_rf <- randomForest(phylum ~., training)
#data2_rf
RF_time <- toc()
#opar <- par(mfrow=c(1,2))
#plot(data2_rf)
#varImpPlot(data2_rf)
#par(opar)

#EVALUATION
yhat_rf <- predict(data2_rf, newdata=validation, type='prob')[,'1']
pred_rf <- prediction(yhat_rf, y_obs)
#performance(pred_rf, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_rf)


#BOOSTING
set.seed(2305)
data2_for_gbm <-
  training %>%
  mutate(phylum=as.numeric(as.character(phylum)))
data2_gbm <- gbm(phylum ~ .,
                 data = data2_for_gbm, distribution="bernoulli",
                 n.trees=131,
                 n.minobsinnode = 1,
                 cv.folds=3,
                 verbose=T)
(best_iter = gbm.perf(data2_gbm, method="cv")) #131

#EVALUATION
yhat_gbm <- predict(data2_gbm, n.trees=best_iter, newdata=validation, type='response')
pred_gbm <- prediction(yhat_gbm, y_obs)
#performance(pred_gbm, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_gbm)


#최종 모형 선택과 테스트세트 오차 계산
data.frame(method=c('lm', 'glmnet', 'rf', 'gbm'),
           auc = c(performance(pred_lm, "auc")@y.values[[1]],
                   performance(pred_glmnet, "auc")@y.values[[1]],
                   performance(pred_rf, "auc")@y.values[[1]],
                   performance(pred_gbm, "auc")@y.values[[1]]),
           bin_dev = c(binomial_deviance(y_obs, yhat_lm),
                       binomial_deviance(y_obs, yhat_glmnet),
                       binomial_deviance(y_obs, yhat_rf),
                       binomial_deviance(y_obs, yhat_gbm)))
# LASSO > boosting > logistic > randomforest

data.frame(time=c(LASSO_time$callback_msg, RF_time$callback_msg))


#ROC curve
perf_lm <- performance(pred_lm, measure = "tpr", x.measure = "fpr")
perf_glmnet <- performance(pred_glmnet, measure = "tpr", x.measure = "fpr")
perf_rf <- performance(pred_rf, measure = "tpr", x.measure = "fpr")
perf_gbm <- performance(pred_gbm, measure = "tpr", x.measure = "fpr")

plot(perf_lm, col='black', main="ROC Curve")
plot(perf_glmnet, add=T, col='blue')
plot(perf_rf, add=T, col='red')
plot(perf_gbm, add=T, col='cyan')
abline(0,1)
legend('bottomright', inset=.1,
       legend=c("GLM", "glmnet", "RF", "GBM"),
       col=c('black', 'blue', 'red', 'cyan'), lty=1, lwd=2)

#VISUALIZATION
#pairs(data.frame(y_obs=y_obs,
#                 yhat_lm=yhat_lm,
#                 yhat_glmnet=c(yhat_glmnet),
#                 yhat_rf=yhat_rf,
#                 yhat_gbm=yhat_gbm),
#      lower.panel=function(x,y){ points(x,y); abline(0, 1, col='red')},
#      upper.panel = panel.cor)

#LASSO test
y_obs_test <- as.numeric(as.character(test$phylum))
yhat_glmnet_test <- predict(data2_cvfit, s="lambda.min", newx=xx[test_idx,],
                            type='response')
yhat_glmnet_test <- yhat_glmnet_test[,1]
pred_glmnet_test <- prediction(yhat_glmnet_test, y_obs_test)
performance(pred_glmnet_test, "auc")@y.values[[1]]
binomial_deviance(y_obs_test, yhat_glmnet_test)
```

##### 2-3. example2: firmicutes & others with collapsed 5mer
```{r example2_data}

# example2 : firmicutes & others with collasped 5 mer
firmicutes_data <- phylum_seq %>% filter(phylum == "firmicutes") %>%
  sample_n(10000)
others_data <- phylum_seq %>% filter(phylum != "firmicutes") %>%
  sample_n(10000)
data <- rbind(firmicutes_data, others_data)


# counting nucleotide pentamer
r16Sstr <- DNAStringSet(data$sequence)
#k1fq = alphabetFrequency(r16Sstr)
#k1fq = k1fq[,1:4]
#k2fq = dinucleotideFrequency(r16Sstr)
#k3fq = trinucleotideFrequency(r16Sstr)
#k4fq = oligonucleotideFrequency(r16Sstr, width=4)
k5fq = oligonucleotideFrequency(r16Sstr, width=5)

#data <- cbind(data[,1], k1fq, k2fq, k3fq, k4fq, k5fq)
data <- cbind(data[,1], k5fq)
colnames(data)[1]<-'phylum'
data$phylum <- factor(ifelse(data$phylum == 'firmicutes', 0, 1))
head(data)

# (optional) Removing highly correlated variables
pairs.panels(data[,c(2:11)], 
             method = "pearson",
             hist.col = "#00AFBB",
             density = TRUE,
             ellipses = TRUE)
library(tictoc)
tic("correlation")
cor_matrix <- cor(data[,c(2:1025)])
cor_time <- toc()
data.frame(time=cor_time$callback_msg)
cor_matrix_rm <- cor_matrix
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
diag(cor_matrix_rm) <- 0

data_new <- data[ , !apply(cor_matrix_rm,
                           2,
                           function(x) any(abs(x) > 0.30))] 

pairs.panels(data_new[,c(1:10)], 
             method = "pearson",
             hist.col = "#00AFBB",
             density = TRUE,
             ellipses = TRUE
)

#data2 <- data
data2 <- cbind(data[,1], data_new)
colnames(data2)[1]<-'phylum'
data.frame(data=dim(data)[2], data2=dim(data2)[2])
```

##### example2 running
```{r running2}
#TRAINING
set.seed(2306)
n <- nrow(data2)
idx <- 1:n
training_idx <- sample(idx, n * .60)
idx <- setdiff(idx, training_idx)
validate_idx <- sample(idx, n * .20)
test_idx <- setdiff(idx, validate_idx)
training <- data2[training_idx,]
validation <- data2[validate_idx,]
test <- data2[test_idx,]


#LOGISTIC REGRESSION
data2_lm_full <- glm(phylum ~ ., data=training, family=binomial)
#summary.glm(data2_lm_full)
#predict(data2_lm_full, newdata = data2[1:5,], type='response')

#MODEL EVALUATION
y_obs <- as.numeric(as.character(validation$phylum))
yhat_lm <- predict(data2_lm_full, newdata = validation, type='response')
pred_lm <- prediction(as.numeric(yhat_lm), as.numeric(y_obs))
#performance(pred_lm, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_lm) 
#plot(performance(pred_lm, measure="tpr", x.measure='fpr'))


#LASSO MODEL
xx <- model.matrix(phylum ~ .-1, data2)
x <- xx[training_idx, ]
y <- as.numeric(as.character(training$phylum))
#glimpse(x)

library(tictoc)
tic("LASSO")
data2_cvfit <- cv.glmnet(x, y, family = "binomial")
#plot(data2_cvfit)
LASSO_time <- toc()

#coef(data2_cvfit, s = c("lambda.1se"))
#coef(data2_cvfit, s = c("lambda.min")) #예측에는 lambda.min 사용

#MODEL EVALUATION
#predict.cv.glmnet(data2_cvfit, s='lambda.min', newx = x[1:5,], type='response')

yhat_glmnet <- predict(data2_cvfit, s="lambda.min", newx=xx[validate_idx,],
                       type='response')
yhat_glmnet <- yhat_glmnet[,1]
pred_glmnet <- prediction(yhat_glmnet, y_obs)
#performance(pred_glmnet, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_glmnet)

#TREE MODEL
data2_tr <- rpart(phylum ~ ., data = training)
#data2_tr

#printcp(data2_tr)
#summary(data2_tr)

opar <- par(mfrow = c(1,1), xpd = NA)
plot(data2_tr)
text(data2_tr, use.n = T)
par(opar)

#EVALUATION
yhat_tr <- predict(data2_tr, validation)
yhat_tr <- yhat_tr[,"1"]
pred_tr <- prediction(yhat_tr, y_obs)
#performance(pred_tr, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_tr)

#RANDOM FOREST
set.seed(2305)
tic("RF")
data2_rf <- randomForest(phylum ~., training)
#data2_rf
RF_time <- toc()
#opar <- par(mfrow=c(1,2))
#plot(data2_rf)
#varImpPlot(data2_rf)
#par(opar)

#EVALUATION
yhat_rf <- predict(data2_rf, newdata=validation, type='prob')[,'1']
pred_rf <- prediction(yhat_rf, y_obs)
#performance(pred_rf, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_rf)


#BOOSTING
set.seed(2305)
data2_for_gbm <-
  training %>%
  mutate(phylum=as.numeric(as.character(phylum)))
data2_gbm <- gbm(phylum ~ .,
                 data = data2_for_gbm, distribution="bernoulli",
                 n.trees=131,
                 n.minobsinnode = 1,
                 cv.folds=3,
                 verbose=T)
(best_iter = gbm.perf(data2_gbm, method="cv")) #131

#EVALUATION
yhat_gbm <- predict(data2_gbm, n.trees=best_iter, newdata=validation, type='response')
pred_gbm <- prediction(yhat_gbm, y_obs)
#performance(pred_gbm, "auc")@y.values[[1]]
#binomial_deviance(y_obs, yhat_gbm)


#최종 모형 선택과 테스트세트 오차 계산
data.frame(method=c('lm', 'glmnet', 'rf', 'gbm'),
           auc = c(performance(pred_lm, "auc")@y.values[[1]],
                   performance(pred_glmnet, "auc")@y.values[[1]],
                   performance(pred_rf, "auc")@y.values[[1]],
                   performance(pred_gbm, "auc")@y.values[[1]]),
           bin_dev = c(binomial_deviance(y_obs, yhat_lm),
                       binomial_deviance(y_obs, yhat_glmnet),
                       binomial_deviance(y_obs, yhat_rf),
                       binomial_deviance(y_obs, yhat_gbm)))
# LASSO > boosting > logistic > randomforest

data.frame(time=c(LASSO_time$callback_msg, RF_time$callback_msg))


#ROC curve
perf_lm <- performance(pred_lm, measure = "tpr", x.measure = "fpr")
perf_glmnet <- performance(pred_glmnet, measure = "tpr", x.measure = "fpr")
perf_rf <- performance(pred_rf, measure = "tpr", x.measure = "fpr")
perf_gbm <- performance(pred_gbm, measure = "tpr", x.measure = "fpr")

plot(perf_lm, col='black', main="ROC Curve")
plot(perf_glmnet, add=T, col='blue')
plot(perf_rf, add=T, col='red')
plot(perf_gbm, add=T, col='cyan')
abline(0,1)
legend('bottomright', inset=.1,
       legend=c("GLM", "glmnet", "RF", "GBM"),
       col=c('black', 'blue', 'red', 'cyan'), lty=1, lwd=2)

#VISUALIZATION
#pairs(data.frame(y_obs=y_obs,
#                 yhat_lm=yhat_lm,
#                 yhat_glmnet=c(yhat_glmnet),
#                 yhat_rf=yhat_rf,
#                 yhat_gbm=yhat_gbm),
#      lower.panel=function(x,y){ points(x,y); abline(0, 1, col='red')},
#      upper.panel = panel.cor)

#LASSO test
y_obs_test <- as.numeric(as.character(test$phylum))
yhat_glmnet_test <- predict(data2_cvfit, s="lambda.min", newx=xx[test_idx,],
                            type='response')
yhat_glmnet_test <- yhat_glmnet_test[,1]
pred_glmnet_test <- prediction(yhat_glmnet_test, y_obs_test)
performance(pred_glmnet_test, "auc")@y.values[[1]]
binomial_deviance(y_obs_test, yhat_glmnet_test)

```

##### 3-1. multivariative by binomial data
```{r example3 data}
# number of phyla which have more than 100, 1000, 10000 freq
phylum_table %>%
  summarize(more_than_100 = n_distinct(phylum_table %>% filter(Freq > 100)),
            more_than_1000 = n_distinct(phylum_table %>% filter(Freq > 1000)),
            more_than_10000 = n_distinct(phylum_table %>% filter(Freq > 10000))
  ) #120, 38, 5


#random 100 sequences of each phylum
filtered_phylum_table <- phylum_table %>% filter(Freq >= 1000)
#filtered_phylum_table2 <- phylum_table %>% filter(Freq <= 100, !is.na(phylum)) # 196 phylum
sampled_phylum_seq <- data.frame(matrix(nrow=0, ncol=2))
colnames(sampled_phylum_seq) <- c('phylum', 'sequence')
for (i in 1:38){
  sampled_phylum_seq <-
    rbind(sampled_phylum_seq, phylum_seq %>%
            filter(phylum == filtered_phylum_table$phylum[i]) %>%
            sample_n(1000))
}
glimpse(sampled_phylum_seq)




phylum_seq_pool <- data.frame(matrix(nrow=0, ncol=2))
for (i in 1:38){
  target <- sampled_phylum_seq %>%
    filter(phylum == filtered_phylum_table$phylum[i])
  learn <- sampled_phylum_seq %>%
    filter(phylum != filtered_phylum_table$phylum[i]) %>%
    sample_n(1000)
  target_learn <- rbind(target, learn)
  phylum_seq_pool <- rbind(phylum_seq_pool, target_learn)
}
dim(phylum_seq_pool)



#sampled_seq_taxa2 <- data.frame(matrix(nrow=0, ncol=2))
#colnames(sampled_seq_taxa2) <- c('phylum', 'sequence')
#for (i in 1:197){
#  sampled_seq_taxa2 <- rbind(sampled_seq_taxa2,
#                             as.data.frame(seq_taxa %>%
#                                             dplyr::select(phylum, sequence) %>%
#                                             filter(phylum == filtered_phylum_table2$phylum[i])))
#}


#sampled_seq_taxa <- rbind(sampled_seq_taxa, sampled_seq_taxa2)

phylum_table %>% filter(Freq >= 1000)

```

##### 3-1 finding
```{r finding}
#example <- phylum_seq %>% filter(phylum == "acidobacteriota") %>%
#  sample_n(1)

example <- readDNAStringSet("/Users/chungminkim/Downloads/sequence-2.fasta")
# e.coli - proteobacteria

example <- data.frame(id=names(example),sequences=as.character(example))
colnames(example) <- c('phylum', 'sequence')

result <- data.frame(matrix(nrow=0, ncol=2))
tic("finding")
for (i in 1:38){
  data <- phylum_seq_pool[(2000*i-1999):(2000*i),]
  #data <- phylum_seq_pool[1:200,]
  data <- rbind(data, example)
  
  r16Sstr <- DNAStringSet(data$sequence)
  k5fq = oligonucleotideFrequency(r16Sstr, width=5)
  
  data <- cbind(data[,1], k5fq)
  colnames(data)[1]<-'phylum'
  #as_tibble(data)
  #data$phylum <- factor(ifelse(data$phylum == phylum_seq_pool[200*i-199,1], 0, 1))
  data[1:1000,1] <- 0
  data[1001:2001,1] <- 1
  cor_matrix <- cor(data[,c(2:1025)])
  cor_matrix_rm <- cor_matrix
  cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
  diag(cor_matrix_rm) <- 0
  
  
  #pairs.panels(data[,c(2:11)], 
  #             method = "pearson",
  #             hist.col = "#00AFBB",
  #             density = TRUE,
  #             ellipses = TRUE)
  
  data_new <- data[ , !apply(cor_matrix_rm,
                             2,
                             function(x) any(abs(x) > 0.30))] 
  
  #pairs.panels(data_new[,c(1:10)], 
  #             method = "pearson",
  #             hist.col = "#00AFBB",
  #             density = TRUE,
  #             ellipses = TRUE
  #)
  
  
  data2 <- cbind(data[,1], data_new)
  colnames(data2)[1]<-'phylum'
  
  example2 <- data2[2001, ]
  data2 <- data2[1:2000, ]
  
  
  #TRAINING
  set.seed(2306)
  n <- nrow(data2)
  idx <- 1:n
  training_idx <- sample(idx, n * .60)
  idx <- setdiff(idx, training_idx)
  validate_idx <- sample(idx, n * .20)
  test_idx <- setdiff(idx, validate_idx)
  training <- data2[training_idx,]
  validation <- data2[validate_idx,]
  test <- data2[test_idx,]
  
  xx <- model.matrix(phylum ~ .-1, data2)
  x <- xx[training_idx, ]
  y <- as.numeric(as.character(training$phylum))
  #glimpse(x)
  
  data2_cvfit <- cv.glmnet(x, y, family = "binomial")
  #plot(data2_cvfit)
  
  #coef(data2_cvfit, s = c("lambda.1se"))
  #coef(data2_cvfit, s = c("lambda.min"))
  
  #MODEL EVALUATION
  
  
  #######
  example2 <- example2[1, 2:(ncol(example2))]
  pre <- as.data.frame(predict.cv.glmnet(data2_cvfit, s='lambda.min',
                                         newx = as.double(example2), type='response'))
  
  result <- rbind(result, pre)
  #######
  
  #yhat_glmnet <- predict(data2_cvfit, s="lambda.min", newx=xx[validate_idx,],
  #                       type='response')
  #yhat_glmnet <- yhat_glmnet[,1]
  #pred_glmnet <- prediction(yhat_glmnet, y_obs)
  #performance(pred_glmnet, "auc")@y.values[[1]]
  #binomial_deviance(y_obs, yhat_glmnet)
}
finding_time <- toc()

#data.frame(time=c(finding_time$callback_msg))

result
phylum_seq_pool[as.vector(which.min(result[,1])) * 2000 - 1999, 1]

rank <- cbind(filtered_phylum_table$phylum, result)
head(arrange(rank, by_group=lambda.min))
```
