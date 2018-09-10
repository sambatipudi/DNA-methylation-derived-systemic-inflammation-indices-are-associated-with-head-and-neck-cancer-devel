## Load relevant data and packages
data <- read.csv("sri_2.csv", stringsAsFactors = F)
data_2 <- read.csv("sri_3.csv", stringsAsFactors = F)
data <- data[,2:ncol(data)]
data$lymph <- (data$Bcell+data$CD4T+data$CD8T+data$NK)
data <- merge(data, data_2, by.x = "IID", by.y = "datafilename", all.x = T)
load("sri_1.RData")
library(reshape2)
library(ggplot2)

setwd("~/mDNLR project/HNSCC_mdNLR")
load("~/mDNLR project/HNSCC_mdNLR/sri_dfs.Rdata")
save(test1, file= "Data_HNSCC_DFS.RData") # The dataset is called test1 
data<- read.csv("data_HNSCC_Ryan.csv", header=T)
colnames(data)
data<-data[,-1]
write.csv(data, file = "data_HNSCC_Ryan.csv")
# Reload the dataset
data<- read.csv("data_HNSCC_Ryan.csv", header=T)
data<-data[,-1]




###################################################################
# Merge test1 and data to get dataset with information on disease free survival
###################################################################
merged<-merge(data, test1, by.x="IID", by.y="datafilename")
write.csv(merged, file = "merged_data_HNSCC_DFS.csv")
# Final_dataset_HNSCC_DFS. csv was made selecting the variables of interest from the merged dataset
# This dataset will be used for all analysis

data<- read.csv("Final_dataset_HNSCC_DFS.csv", header=T)
out <- strsplit(as.character(data$IID),'_') 
do.call(rbind, out)
df<-data.frame(data$IID, do.call(rbind, out))
merge<-merge(data, df, by.x="IID", by.y="data.IID")
colnames(merge)
colnames(merge)[25] <- "Sentrix_ID"
colnames(merge)[26] <- "Sentrix_Pos"
merge<-merge[,-1]
library(dplyr)
new_df <- merge %>% select(Sentrix_ID, Sentrix_Pos, everything())

save(new_df, file= "Data_Final_dataset_HNSCC_DFS.RData")


load("~/mDNLR project/HNSCC_mdNLR/Data_Final_dataset_HNSCC_DFS.RData")

data<-new_df

data$DFS_1 <- ifelse(data$DFS=="Curative", 0,1)

save(data_mdNLR, file= "Data_Final_dataset_HNSCC_DFS_mdNLR.RData")

save(data_mdLMR, file= "Data_Final_dataset_HNSCC_DFS_mdLMR.RData")

write.csv(data, file = "Data_HNSCC_DFS.csv")

# Replacing DFS values with 0 and 1 for DFS analysis
data$DFS <- ifelse(data$DFS=="Curative", 0,1)



## Boxplots in ggplot2 (Early, late stage, HPV+, HPVNeg, Dead , Alive etc.)
box_stage <- data[,c(17,6, 7, 18)]
colnames(box_stage) <- c("Stage","Monocytes", "Granulocytes", "Lymphocytes")

box_hpv <- data[,c(16, 6, 7, 18)]
colnames(box_hpv) <- c("HPV16E6","Monocytes", "Granulocytes", "Lymphocytes")

box_doa <- data[,c(19, 6, 7, 18)]
colnames(box_doa) <- c("Survival_status", "Monocytes", "Granulocytes", "Lymphocytes")

box_smoke <- data[,c(13, 6, 7, 18)]
colnames(box_smoke) <- c("Smoking","Monocytes", "Granulocytes", "Lymphocytes")
# There are NAs in Smoking staus column (Remove them from dataframe)
table(is.na(box_smoke$Smoking))
box_smoke<-box_smoke[complete.cases(box_smoke), ]


box_drink <- data[,c(14, 6, 7, 18)]
colnames(box_drink) <- c("Drinking","Monocytes", "Granulocytes", "Lymphocytes")


box_gender <- data[,c(12, 6, 7, 18)]
colnames(box_gender) <- c("Gender","Monocytes", "Granulocytes", "Lymphocytes")


box_age <- data[,c(11, 6, 7, 18)]
colnames(box_age) <- c("DichotomisedAge", "Monocytes", "Granulocytes", "Lymphocytes")

# plot the proportions by group of interest
colnames(box_stage)
box_stage$Stage[box_stage$Stage == 0] <- "Low"
box_stage$Stage[box_stage$Stage == 1] <- "High"

stage.m <- melt(box_stage)
tiff("boxplot_stage.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Stage)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("mdNLR")+
  ggtitle("")
dev.off()



box_hpv$HPV16E6[box_hpv$HPV16E6 == 0] <- "Negative"
box_hpv$HPV16E6[box_hpv$HPV16E6 == 1] <- "Positive"

stage.m <- melt(box_hpv)
tiff("boxplot_HPV.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = HPV16E6)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("mdNLR")+
  ggtitle("")
dev.off()



box_doa$Survival_status[box_doa$Survival_status == 0] <- "Alive"
box_doa$Survival_status[box_doa$Survival_status == 1] <- "Dead"

stage.m <- melt(box_doa)
tiff("boxplot_survial.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Survival_status)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("mdNLR")+
  ggtitle("")
dev.off()


box_smoke$Smoking[box_smoke$Smoking == 0] <- "Never"
box_smoke$Smoking[box_smoke$Smoking == 1]  <- "Ever"

stage.m <- melt(box_smoke)
tiff("boxplot_smoking.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Smoking)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("mdNLR")+
  ggtitle("")
dev.off()


box_drink$Drinking[box_drink$Drinking == 0] <- "Never"
box_drink$Drinking[box_drink$Drinking == 1]  <- "Ever"

stage.m <- melt(box_drink)
tiff("boxplot_drinking.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Drinking)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("Estimated proportion in whole blood")+
  ggtitle("")
dev.off()


box_gender$Gender[box_gender$Gender == 0] <- "Females"
box_gender$Gender[box_gender$Gender == 1]  <- "Males"

stage.m <- melt(box_gender)
tiff("boxplot_gender.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Gender)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("Estimated proportion in whole blood")+
  ggtitle("")
dev.off()



box_age$DichotomisedAge[box_age$DichotomisedAge == 0] <- "<60 Years"
box_age$DichotomisedAge[box_age$DichotomisedAge == 1]  <- ">=60 Years"

stage.m <- melt(box_age)
tiff("boxplot_age.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = DichotomisedAge)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("Estimated proportion in whole blood")+
  ggtitle("")
dev.off()


# Neutrophil to Lymphocyte Ratio (NLR) analysis
# 1. Age (No)
lm<-lm(mdNLR~Age, data=data) 
summary(lm)

# 2. Gender (No)
lm<-lm(mdNLR~gender, data=data)
summary(lm)


# 3. Drinking status (continous) (NO)
lm<-lm(mdNLR~en_drink, data=data)
summary(lm)

# 4. Smoking status (Ever/Never) (Yes, marginal)P= 0.057 .
lm<-lm(mdNLR~en_smoke+as.factor(Sentrix_ID), data=data)
summary(lm)

cbind(coef(lm), confint(lm))

# 5. Tumour grade (Low/High) (Yes, marginal, P=0.0572)

lm<-lm(mdNLR~grade+as.factor(Sentrix_ID), data=data)
summary(lm)

cbind(coef(lm), confint(lm))


lm<-lm(mdNLR~grade+gender, data=data) #  P=0.0634  
summary(lm)

lm<-lm(mdNLR~grade+gender+en_smoke+en_drink, data=data) #P=0.0520
summary(lm)

lm<-lm(mdNLR~grade+gender+en_smoke+en_drink+Age, data=data)#P=0.048 *
summary(lm)

lm<-lm(mdNLR~grade+gender+en_smoke+en_drink+Age+HPV16E6_pos2+as.factor(Sentrix_ID), data=data)# P=0.041 *
summary(lm)
cbind(coef(lm), confint(lm))



# 6. HPV status (No)
lm<-lm(mdNLR~HPV16E6_pos2+as.factor(Sentrix_ID), data=data)
summary(lm)


# 7. Sentrix_ID of samples (Associated with some Sentrix IDs)

lm<-lm(mdNLR~as.factor(Sentrix_ID), data=data)
summary(lm)

####################################################################

# Uncoditional logistic regression model for overll and disease free survival

####################################################################
# Disease Free Survival
# Univariate analysis # 12 samples with NA in DFS removed
fit1<- glm(DFS ~ mdNLR, family=binomial(link='logit'),data=data)
summary(fit1)

exp(cbind(OR=coef(fit1), confint(fit1)))


#library(pROC)
#a=rbinom(100, 1, 0.25)
#b=runif(100)
#c=rnorm(100)

#fit1=glm(a~b+c, family='binomial')
#fit2=glm(a~c, family='binomial')

#preds=predict(fit1)
#roc1=roc(data ~ preds)

#preds2=predict(fit2)
#roc2=roc(a ~ preds2)

#plot(roc1, print.auc = TRUE)
#plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)



#### Plotting AUC with data on mdNLR and mdLMR
# Using the whole dataset for mdNLR: Disease free survival

#fit2<- glm(DFS ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
#summary(fit2)

fit2<- glm(DFS ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(fit2)

fit3<- glm(DFS ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
summary(fit3)

exp(cbind(OR=coef(fit3), confint(fit3)))


merged<-merge(data, test1, by.x="IID", by.y="datafilename")



# Using the whole dataset for mdNLR: Overall survival
# Univariate model

fit1<- glm(ysnPatientDead ~ mdNLR, family=binomial(link='logit'),data=data)
summary(fit1)
exp(cbind(OR=coef(fit1), confint(fit1)))

#fit2<- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(fit2)

fit3<- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(fit3)

fit4<- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
summary(fit4)
options(scipen=999)
options(digits=2)
exp(cbind(OR=coef(fit4), confint(fit4)))



preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

plot(roc1,col='green', print.auc = TRUE)
plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = TRUE, print.auc.y = .3)

# Using the whole dataset for mdLMR: Disease free survival

lm<-lm(mdLMR~as.factor(Sentrix_ID), data=data)
summary(lm)


# Disease Free Survival
# Univariate analysis # 

fit1<- glm(DFS ~ mdLMR, family=binomial(link='logit'),data=data)
summary(fit1)
exp(cbind(OR=coef(fit1), confint(fit1)))


#fit2<- glm(DFS ~ Age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)


# Multivariate analysis (without Sentrx_ID)

fit3<- glm(DFS ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(fit3)

fit4<- glm(DFS~  mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
summary(fit4)
exp(cbind(OR=coef(fit4), confint(fit4)))


# Using the whole dataset for mdLMR: Overall survival

fit1<- glm(ysnPatientDead ~ mdLMR, family=binomial(link='logit'),data=data)
summary(fit1)
exp(cbind(OR=coef(fit1), confint(fit1)))

fit2<- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(fit2)


fit3<- glm(ysnPatientDead ~  mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
summary(fit3)
exp(cbind(OR=coef(fit3), confint(fit3)))

##########################################################
# AUC curve for prediction of survival for mdNLR and mdLMR
##########################################################
#DFS
fit1<- glm(DFS ~ mdNLR, family=binomial(link='logit'),data=data)
fit2<- glm(DFS ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
fit3<- glm(DFS ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

tiff("ROC_modelcomparison_DFS_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdNLR only"=roc1, "Covariates only"=roc2, "mdNLR+Covariates"=roc3),aes="colour", size=2)
dev.off()

tiff("ROC_modelcomparison_DFS_NLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2, print.auc.x = .7, print.auc.y = .3, xlim=c(1.0, 0.0))
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.x = .7,print.auc.y = .2, xlim=c(1.0, 0.0))
plot(roc3, add=TRUE, col='blue', print.auc = T, print.auc.x = .7,print.auc.y = .1, xlim=c(1.0, 0.0))
dev.off()

#Overall survival
fit1<- glm(ysnPatientDead ~ mdNLR, family=binomial(link='logit'),data=data)
fit2<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
fit3<- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

tiff("ROC_modelcomparison_OS_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdNLR only"=roc1, "Covariates only"=roc2, "mdNLR+Covariates"=roc3),aes="colour", size=2)
dev.off()


tiff("ROC_modelcomparison_OS_mdNLR_II.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2,  print.auc.x = .7, print.auc.y = .3)
plot(roc2, add=TRUE, col='red', print.auc = T,  print.auc.x = .7, print.auc.y = .2)
plot(roc3, add=TRUE, col='blue', print.auc = T,  print.auc.x = .7, print.auc.y = .1)
dev.off()



### #Overall survival adjusted for Sentrix_ID
fit1<- glm(ysnPatientDead ~ mdNLR +as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
fit2<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
fit3<- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)


tiff("ROC_modelcomparison_OS_mdNLR_adjSentrix.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdNLR only"=roc1, "Covariates only"=roc2, "mdNLR+Covariates"=roc3),aes="colour", size=2)
dev.off()

tiff("ROC_modelcomparison_OS_mdNLR_II_adjSentrix.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = T, print.auc.y = .3)
dev.off()



#DFS
fit1<- glm(DFS ~ mdLMR, family=binomial(link='logit'),data=data)
fit2<- glm(DFS ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
fit3<- glm(DFS ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

tiff("ROC_modelcomparison_DFS_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdLMR only"=roc1, "Covariates only"=roc2, "mdLMR+Covariates"=roc3),aes="colour", size=2)
dev.off()

tiff("ROC_modelcomparison_DFS_mdLMR_II.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = T, print.auc.y = .3)
dev.off()

#Overall survival
fit1<- glm(ysnPatientDead ~ mdLMR, family=binomial(link='logit'),data=data)
fit2<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
fit3<- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)




tiff("ROC_modelcomparison_OS_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdLMR only"=roc1, "Covariates only"=roc2, "mdLMR+Covariates"=roc3),aes="colour", size=2)
dev.off()

tiff("ROC_modelcomparison_OS_mdLMR_II.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = T, print.auc.y = .3)
dev.off()

### #Overall survival adjusted for Sentrix_ID
fit1<- glm(ysnPatientDead ~ mdLMR +as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
fit2<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
fit3<- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)




tiff("ROC_modelcomparison_OS_mdLMR_adjSentrix.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdLMR only"=roc1, "Covariates only"=roc2, "mdLMR+Covariates"=roc3),aes="colour", size=2)
dev.off()

tiff("ROC_modelcomparison_OS_mdLMR_II_adjSentrix.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = T, print.auc.y = .3)
dev.off()



### Comapring mdNLR and mdLMR in DFS and OS


### #Disease free survival adjusted for Sentrix_ID
fit1<- glm(DFS ~ mdNLR +as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
fit2<- glm(DFS ~ mdLMR+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)




tiff("ROC_mdNLR_mdLMR_modelcomparison_DFS_adjSentrix.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdNLR"=roc1, "mdLMR"=roc2),aes="colour", size=2)
dev.off()

tiff("ROC_mdNLR_mdLMR_modelcomparison_DFS_adjSentrix_II.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2, print.auc.y = .3)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .2)
dev.off()

### #Overall survival adjusted for Sentrix_ID
fit1<- glm(ysnPatientDead ~ mdNLR +as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
fit2<- glm(ysnPatientDead ~ mdLMR+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)




tiff("ROC_mdNLR_mdLMR_modelcomparison_OS_adjSentrix.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("mdNLR"=roc1, "mdLMR"=roc2),aes="colour", size=2)
dev.off()

tiff("ROC_mdNLR_mdLMR_modelcomparison_OS_adjSentrix_II.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='black', print.auc = T, lwd=2, print.auc.y = .3)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .2)
dev.off()






########################################################################################################



preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

preds4=predict(fit4)
roc4=roc(fit4$y , fit4$fitted.values,ci=T)



plot(roc1,col='green', print.auc = F, lwd=2)
plot(roc2, add=TRUE, col='red', print.auc = T, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = T, print.auc.y = .3)
plot(roc4, add=TRUE, col='black', print.auc = T, print.auc.y = .2)
# Comparing mdNLR and mdLMR for predicting death
fit1<- glm(Status ~ mdNLR, family=binomial(link='logit'),data=data)
fit2<- glm(Status ~ mdLMR, family=binomial(link='logit'),data=data)
preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

plot(roc1,col='green', print.auc = TRUE)
plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)






# Switch off the scientific calculator
options(scipen=999)
options(digits=2)
# Unconditional logistic regression
# Univariate model
model1<- glm(ysnPatientDead ~ mdNLR,family=binomial(link='logit'),data=data)
summary(model1)

# Multivariate model (gender)
model2 <- glm(ysnPatientDead ~ mdNLR+as.factor(gender), family=binomial(link='logit'),data=data)
summary(model2)

# Multivariate model (Age, gender)
model3 <- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender), family=binomial(link='logit'),data=data)
summary(model3)

# Multivariate model (Age, gender, Drinking status, Smoking status, Tumour grade, HPV16E6
# seropositivity)
model4 <- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(en_drink)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(model4)
exp(cbind(coef(model4), confint(model4))) 

# Multivariate model (Age, gender, Smoking status, Tumour grade, HPV16E6
# seropositivity, Sentrix_ID) 
model5 <- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
summary(model5)
exp(cbind(coef(model5), confint(model5)))


# Multivariate model (Smoking status, Tumour grade)- Based on the lm association between mdNLR and variables
#model6 <- glm(ysnPatientDead ~ mdNLR+as.factor(en_smoke)+as.factor(grade), family=binomial(link='logit'),data=data)
#summary(model6)
#exp(cbind(coef(model6), confint(model6))) 



# Lymphocyte to Monocyte Ratio (LMR) analysis
# 1. Age (Yes)
lm<-lm(mdLMR~Age+as.factor(Sentrix_ID), data=data) #P= 0.000033 ***
summary(lm)

# 2. Gender (Yes)
lm<-lm(mdLMR~gender+as.factor(Sentrix_ID), data=data)#P=0.00022 ***
summary(lm)


# 3. Drinking status (continous) (No)
lm<-lm(mdLMR~en_drink+as.factor(Sentrix_ID), data=data)#P=0.45
summary(lm)

# 4. Smoking status (Ever/Never) (Yes)
lm<-lm(mdLMR~en_smoke+as.factor(Sentrix_ID), data=data)#P=0.032 * 
summary(lm)

# 5. Tumour grade (Low/High) No
lm<-lm(mdLMR~grade+as.factor(Sentrix_ID), data=data)#P=0.12 
summary(lm)

lm<-lm(mdLMR~grade+gender, data=data)
summary(lm)

lm<-lm(mdLMR~grade+gender+en_smoke+en_drink, data=data)
summary(lm)

lm<-lm(mdLMR~grade+gender+en_smoke+en_drink+Age, data=data)
summary(lm)

lm<-lm(mdLMR~as.factor(grade)+as.factor(gender)+as.factor(en_smoke)+en_drink+Age+as.factor(HPV16E6_pos2), data=data)
summary(lm)

lm<-lm(mdLMR~as.factor(grade)+as.factor(gender)+as.factor(en_smoke)+en_drink+Age+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data=data)
summary(lm)
cbind(coef(lm), confint(lm)) 



lm<-lm(mdNLR~as.factor(grade)+as.factor(gender)+as.factor(en_smoke)+en_drink+Age+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data=data)
summary(lm)
cbind(coef(lm), confint(lm)) 

# 6. HPV status (Yes)
lm<-lm(mdLMR~HPV16E6_pos2+as.factor(Sentrix_ID), data=data)
summary(lm)

cbind(coef(lm), confint(lm)) 

# 7. Barcode of samples (Did not work)

lm<-lm(mdLMR~IID, data=data)
summary(lm)


# Unconditional logistic regression
# Univariate model
model1<- glm(ysnPatientDead ~ mdLMR,family=binomial(link='logit'),data=data)
summary(model1)

# Multivariate model (gender)
model2 <- glm(ysnPatientDead ~ mdLMR+as.factor(gender), family=binomial(link='logit'),data=data)
summary(model2)

# Multivariate model (Age, gender)
model3 <- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender), family=binomial(link='logit'),data=data)
summary(model3)

# Multivariate model (Age, gender, Smoking status, Tumour grade, HPV16E6
# seropositivity)
model4 <- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), family=binomial(link='logit'),data=data)
summary(model4)
exp(cbind(coef(model4), confint(model4))) 


# Multivariate model (Age, gender, Drinking status, Smoking status, Tumour grade, HPV16E6
# seropositivity, Sentrix_ID)
model5 <- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), family=binomial(link='logit'),data=data)
summary(model5)
exp(cbind(coef(model5), confint(model5))) 
##################################################
# Survival analysis AUC analysis using ROCR
##################################################
library(ROCR)
colnames(data)
data$Y = as.numeric(data$ysnPatientDead == "1")
pred1 <- prediction(data$mdNLR, data$Y)
perf1 <- performance(pred1,"tpr","fpr")
tiff("ROC_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(perf1)
abline(a=0, b= 1)
dev.off()

opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))

#              [,1]
#sensitivity 0.8494624
#specificity 0.8504673
#cutoff      2.7214037


data$Y = as.numeric(data$ysnPatientDead == "0")
pred2 <- prediction(data$mdLMR, data$Y)
perf2 <- performance(pred2,"tpr","fpr")
tiff("ROC_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(perf2)
abline(a=0, b= 1)
dev.off()

opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))
#[,1]
#sensitivity 0.8494624
#specificity 0.8504673
#cutoff      4.4528157

acc.perf = performance(pred, measure = "acc")
tiff("Accuracy_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(acc.perf)
dev.off()


#############################################   
# PLotting mdNLR and mdLMR in the same plot
#############################################

tiff("AUC_mdNLR_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot( perf1, colorize = TRUE)
plot(perf2, add = TRUE, colorize = TRUE, xlab="", ylab="")
abline(a=0, b= 1)
dev.off()


range(data$age)

# Load mdNLR data for HNSCC

load("//ads.bris.ac.uk/filestore/MyFiles/Staff16/sa16666/Documents/mDNLR project/HNSCC_mdNLR/boxplots.RData")
dim(data)
colnames(data)

# Correlation between mdNLR and mdLMR in samples
summary(SABRE_merge_base$bl_crpv1)
cor(data$mdNLR, data$mdLMR, use = "complete.obs")
cor.test(data$mdNLR, data$mdLMR, use = "complete.obs") # p-value < 2.2e-16

library("ggpubr")
tiff("Corr_mdNLR_mdLMR.tif", res=300, compression = "lzw", height=8, width=10, units="in")
ggscatter(data, x = "mdLMR" , y = "mdNLR", add = "reg.line", conf.int = TRUE,ylim = c(0,10), 
          cor.coef = TRUE, cor.method = "spearman",cor.coef.coord = c(8,10),cor.coef.size = 5,
          xlab = "mdLMR", ylab = "mdNLR")
dev.off()


#####################################

# Read the datasent by Ryan 070917
# STATA file

#####################################
setwd("~/mDNLR project/HNSCC_mdNLR")
library(foreign)
library(readstata13)
mydata <- read.dta("HN003_H&N dataset_v2.1_stata12.dta")

merged.data <- merge(data, mydata, by.x="IID", by.y="datafilename")

save(merged.data, file="HNSCC5000_mdNLR_090917.RData")

write.csv(merged.data, file="HNSCC5000_mdNLR_090917.csv")

# The merged data was modified to have only varaibles of interest rest was removed from the df
# Load the moidifed data farme
data<-read.csv("HNSCC5000_mdNLR_090917_mod.csv")

head(data)


data("lung")
head(lung)
res.cox <- coxph(Surv(time_months, Status) ~ mdNLR, data = data)
res.cox

res.cox <- coxph(Surv(time_months, Status) ~ mdNLR+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(en_drink)+as.factor(HPV16E6_pos2.x), data = data)
res.cox


res.cox <- coxph(Surv(time_months, Status) ~ mdNLR, data = data)
res.cox


res.cox <- coxph(Surv(time_months, Status) ~ mdLMR+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(en_drink)+as.factor(HPV16E6_pos2.x), data = data)
res.cox


res.cox <- coxph(Surv(time_months, Status) ~ data$mdNLR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(en_drink)+as.factor(HPV16E6_pos2.x), data = data)
res.cox

res.cox <- coxph(Surv(time_months, Status) ~ data$mdNLR_ab_bel, data = data)
res.cox

res.cox <- coxph(Surv(time_months, Status) ~ 1, data = data)
res.cox

cox_fit <- survfit(res.cox)

summary(cox_fit)
plot(cox_fit)


res.cox <- coxph(Surv(time_months, Status) ~ mdNLR_ab_bel , data = data)
res.cox

res.cox <- survfit(Surv(time_months, Status) ~ mdNLR_ab_bel , data = data)
res.cox


res.cox <- survfit(Surv(time_months, Status) ~ mdNLR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data = data)

summary(res.cox)

res.cox <- coxph(Surv(time_months, Status) ~ mdNLR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(en_drink)+as.factor(HPV16E6_pos2.x), data = data)

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))

# testing the validity of CPH model
test.ph <- cox.zph(res.cox)

test.ph

ggcoxzph(test.ph)


### Removing outliers for mdNLR dataset
newdata <- data[ which(data$mdNLR<= 1.900941),]
                         
res.cox <- coxph(Surv(time_months, Status) ~ mdNLR+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data = newdata)


test.ph <- cox.zph(res.cox)

test.ph

exp(cbind(HR=coef(res.cox), confint(res.cox)))

newdata$mdNLR_ab_bel <- ifelse(newdata$mdNLR<= 1.34681, "Below Median",
                            ifelse(newdata$mdNLR> 1.34681, "Above Median",NA  ))

res.cox <- coxph(Surv(time_months, Status) ~ mdNLR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data = newdata)

exp(cbind(HR=coef(res.cox), confint(res.cox)))

test.ph <- cox.zph(res.cox)

test.ph


### Removing outliers for mdLMR dataset
IQR(data$mdLMR) #IQR=1.765176
1.5*1.765176 # Values above this outliers= 2.647764

newdata_mdLMR <- data[ which(data$mdLMR<= 2.647764),]

res.cox <- coxph(Surv(time_months, Status) ~ mdLMR+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data = newdata)


test.ph <- cox.zph(res.cox)

test.ph

range (data$mdNLR)


sd(data$mdNLR)

range (data$mdLMR)

sd(data$mdLMR)

plot(data$mdNLR)
library(outliers)

#summary(scores(data$mdNLR, type="chisq", prob=0.95))

summary(scores(data$mdNLR, type="z", prob=0.95))

#summary(scores(data$mdNLR, type="t", prob=0.95))

summ=data.frame(scores(data$mdNLR, type="z", prob=0.95))
summ$ID<-row.names(summ)
outlier_ID<-summ[summ$scores.data.mdNLR..type....z...prob...0.95.=="TRUE",]
ID<-as.character(outlier_ID$ID)
as.character(ID)
data_mdNLR<-data[!rownames(data)%in%ID,]


# After removing outliers
library(survival)

# Univariate model

res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR, data = data_mdNLR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data_mdNLR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Including interaction between smoking and mdNLR ( I don't understand why?)
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR*as.factor(en_smoke)+Age+as.factor(gender)+as.factor(grade)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data_mdNLR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))


# Before removing outliers

# Univariate model (CPH assumptions not violated)
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR, data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

# Disease free survival analysis (DFS)
res.cox <- coxph(Surv(TTD_months, DFS) ~ mdNLR, data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

#Multivariate model


res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data)
test.ph <- cox.zph(res.cox)

test.ph
res.cox

#Multivariate model DFS (CPH was not violated, No need to stratify)


res.cox <- coxph(Surv(TTD_months, DFS) ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data)
test.ph <- cox.zph(res.cox)

test.ph

res.cox

#Multivariate model DFS (CPH was violated, stratify by age)
# mdNLR above and below median


res.cox <- coxph(Surv(TTD_months, DFS) ~ as.factor(mdNLR_ab_bel)+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data)
test.ph <- cox.zph(res.cox)

test.ph

res.cox

# Stratify by Age

res.cox <- coxph(Surv(TTD_months, DFS) ~ as.factor(mdNLR_ab_bel)+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID)+strata(Age), data = data)
test.ph <- cox.zph(res.cox)

test.ph
res.cox



# CPH assumption violation for grade and HPV status
# Stratifying the model by grade

res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID)+strata(as.factor(grade)), data = data)

test.ph <- cox.zph(res.cox)

test.ph

# CPH assumption still violated for HPV status


# Stratifying by HPV status


res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(grade)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2)), data = data)

test.ph <- cox.zph(res.cox)

test.ph

# CPH assumption still violated for grade

# Stratifying by HPV status and grade (https://stat.ethz.ch/education/semesters/ss2011/seminar/contents/presentation_5.pdf)


res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))


res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(en_smoke)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data)

res.cox



#### Testing the same model for outlier removed samples
# Multivariate Model
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdNLR+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data_mdNLR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))



#### Outliers for mdLMR 

# Before removing outliers

# Univariate model 
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdLMR, data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox


# Univariate model DFS
res.cox <- coxph(Surv(TTD_months, DFS) ~ mdLMR, data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox




res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))


# Multivariate model DFS (CPH assumptions not violated)
res.cox <- coxph(Surv(TTD_months, DFS) ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))

ggcoxadjustedcurves(res.cox, data = data)


# Multivariate model DFS mdLMR above and below median (Age violated CPH)
res.cox <- coxph(Surv(TTD_months, DFS) ~ as.factor(mdLMR_ab_bel)+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))
fit <- survfit(Surv(TTD_months, DFS) ~ as.factor(mdLMR_ab_bel)+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID), data = data)

# Visualize with survminer
plot(fit)

library(survMisc)
survMisc:::autoplot.survfit(fit)
# Multivariate model DFS mdLMR above and below median (Age violated CPH)
res.cox <- coxph(Surv(TTD_months, DFS) ~ as.factor(mdLMR_ab_bel)+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+as.factor(Sentrix_ID)+strata(Age), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))


# Grade and HPV status are not following the CPH assumption so stratfy the model by grade and HPV status similar to NLR

# Full dataset
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdLMR+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox

exp(cbind(HR=coef(res.cox), confint(res.cox)))




# Removing the Sentrix_ID from the model

res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdLMR+Age+as.factor(gender)+as.factor(en_smoke)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox



##
#####################################################################
summ=data.frame(scores(data$mdLMR, type="z", prob=0.95))
summ$ID<-row.names(summ)
outlier_ID<-summ[summ$scores.data.mdLMR..type....z...prob...0.95.=="TRUE",]
ID<-as.character(outlier_ID$ID)
as.character(ID)
data_mdLMR<-data[!rownames(data)%in%ID,]

# 41 samples removed
####################################################################

# After removing outliers
# Univariate model
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdLMR, data = data_mdLMR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate model
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ mdLMR+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data_mdLMR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))



# add drinking ( Does not change much, so do not include drinking in the model)
res.cox <- coxph(Surv(time_months, Status) ~ mdLMR+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(en_drink)+as.factor(HPV16E6_pos2.x), data = data_mdLMR)

test.ph <- cox.zph(res.cox)

test.ph

res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))


# 

median(data$mdLMR) # 3.2

data$mdLMR_ab_bel <- ifelse(data$mdLMR<= 3.2, "Below Median",
                            ifelse(data$mdLMR> 3.2, "Above Median",NA  ))

median(data$mdNLR) # 1.9

data$mdNLR_ab_bel <- ifelse(data$mdNLR<= 1.9, "Below Median",
                            ifelse(data$mdNLR> 1.9, "Above Median",NA  ))

# Above and below median mdNLR
# Multivariate model
res.cox <- coxph(Surv(TTD_months, ysnPatientDead) ~ as.factor(mdNLR_ab_bel)+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data)

test.ph <- cox.zph(res.cox)

test.ph

res.cox




p <- plot(survfit(res.cox), ylim=c(.45, 1), xlab="Days since treatment start",
          mark.time=T, ylab="Probability of a successful outcome", col=c(1, 2),
          main="Cox proportional hazard model by initial fluoroquinolone
          resistance")

plot(p)



exp(cbind(HR=coef(res.cox), confint(res.cox)))

# survfit <- survfit(Surv(TTD_months, ysnPatientDead) ~ mdNLR_ab_bel+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(Sentrix_ID)+strata(as.factor(HPV16E6_pos2),as.factor(grade)), data = data)

ggsurvplot(survfit, data = data, risk.table = TRUE)


library(survminer)
fit <- survfit(res.cox, newdata = new_df)
ggsurvplot(survfit, conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, surv.median.line = "hv")







survdiff(Surv(time_months, Status) ~ mdLMR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data=data_mdLMR, rho= 0)

survdiff(Surv(time_months, Status) ~ mdLMR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data=data_mdLMR, rho= 1)

survdiff(Surv(time_months, Status) ~ mdLMR_ab_bel+age+as.factor(gender)+as.factor(Stage)+as.factor(en_smoke)+as.factor(HPV16E6_pos2.x), data=data_mdLMR, rho= 1.5)


as.factor(HPV16E6_pos2.x)

# Create the new data  
new_df <- with(data_mdLMR,
               data.frame(mdLMR_ab_bel = c("Above Median", "Below Median"), 
                          age = rep(mean(age, na.rm = TRUE), 2),
                          gender = c("Female", "Male"),
                          Stage = c(0, 1),
                          smoke = c("Ever", "Never"),
                          HPV16E6_pos2.x = c(0,1),
                          )
)
new_df


fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
           ggtheme = theme_minimal())



res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

ggsurvplot(
  survfit,                     # survfit object with calculated statistics.
  data = data,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,70),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)


ggcoxadjustedcurves(res.cox, data = data, variable= data[,"mdNLR_ab_bel"])


ggcoxadjustedcurves(res.cox, data = data,
                    variable  = data [, "mdNLR_ab_bel"],   # Variable of interest
                    legend.title = "mdNLR",        # Change legend title
                     confint=TRUE,   
                     palette = "npg",             # nature publishing group color palettes
                    curv.size = 2                # Change line size
)



ggforest(res.cox, data = survfit, main = "Hazard ratio") 


tiff("survivalplot_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(survfit, conf.int=TRUE,col=c('blue', 'red'), xlab="Months", ylab="Proportion surviving")
legend('bottomleft', c("mdLMR ??? median", "mdLMR > median"), col=c('blue','red'), ncol = 1,
       cex = 0.75, lty=1)
dev.off()
#############################################

# Merging the sva and dataframe 

#############################################
library(tidyr)
unite<-unite(data,basename,Sentrix_ID,Sentrix_Pos, sep = "_")
sva<-read.table("sva_pheno_final.txt", header=T)
merge_sva_data<-merge(unite, sva, by.x = "basename", by.y = "me")
myeloid_CpGs<-read.table ("myeloid.txt", header=T)
colnames(myeloid_CpGs) <- gsub("^X", "",  colnames(myeloid_CpGs))
tmyeloid_CpGs<-data.frame(t(myeloid_CpGs))
tmyeloid_CpGs$basename<-rownames(tmyeloid_CpGs)
merge_tmyeloid_CpGs<-merge(tmyeloid_CpGs, unite, by = "basename")

res.by <- by( merge_tmyeloid_CpGs$cg25938803, merge_tmyeloid_CpGs$ysnPatientDead, mean)
res.by
options(scipen=999)
cg25938803_good<-merge_tmyeloid_CpGs$cg25938803[merge_tmyeloid_CpGs$ysnPatientDead=="0"]
cg25938803_poor<-merge_tmyeloid_CpGs$cg25938803[merge_tmyeloid_CpGs$ysnPatientDead=="1"]
wilcox.test(cg25938803_good,cg25938803_poor)#  0.00009772

res.by <- by( merge_tmyeloid_CpGs$cg10456459, merge_tmyeloid_CpGs$ysnPatientDead, mean)
res.by
cg10456459_good<-merge_tmyeloid_CpGs$cg10456459[merge_tmyeloid_CpGs$ysnPatientDead=="0"]
cg10456459_poor<-merge_tmyeloid_CpGs$cg10456459[merge_tmyeloid_CpGs$ysnPatientDead=="1"]
wilcox.test(cg10456459_good,cg10456459_poor)#0.00006613


res.by <- by( merge_tmyeloid_CpGs$cg01591037, merge_tmyeloid_CpGs$ysnPatientDead, mean)
res.by

cg01591037_good<-merge_tmyeloid_CpGs$cg01591037[merge_tmyeloid_CpGs$ysnPatientDead=="0"]
cg01591037_poor<-merge_tmyeloid_CpGs$cg01591037[merge_tmyeloid_CpGs$ysnPatientDead=="1"]
wilcox.test(cg01591037_good,cg01591037_poor)#0.002576

res.by <- by( merge_tmyeloid_CpGs$cg03621504, merge_tmyeloid_CpGs$ysnPatientDead, mean)
res.by

cg03621504_good<-merge_tmyeloid_CpGs$cg03621504[merge_tmyeloid_CpGs$ysnPatientDead=="0"]
cg03621504_poor<-merge_tmyeloid_CpGs$cg03621504[merge_tmyeloid_CpGs$ysnPatientDead=="1"]
wilcox.test(cg03621504_good,cg03621504_poor)#0.01142

res.by <- by( merge_tmyeloid_CpGs$cg00901982, merge_tmyeloid_CpGs$ysnPatientDead, mean)
res.by
cg00901982_good<-merge_tmyeloid_CpGs$cg00901982[merge_tmyeloid_CpGs$ysnPatientDead=="0"]
cg00901982_poor<-merge_tmyeloid_CpGs$cg00901982[merge_tmyeloid_CpGs$ysnPatientDead=="1"]
wilcox.test(cg00901982_good,cg00901982_poor)# 0.0004954

# Plot the CpGs
test<-merge_tmyeloid_CpGs[,c(24,2:6)]
test$ysnPatientDead<-as.factor(test$ysnPatientDead)
head(test)

test$ysnPatientDead<-ifelse(test$ysnPatientDead==0,"Alive","Dead")

library(reshape2)
test.m <- melt(test, id= "ysnPatientDead")
library(ggplot2)
tiff("boxplot_myeloid_CpGs_HNSCC.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = ysnPatientDead)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("mdNLR associted CpGs") +
  ylab("Methylation")+
  ggtitle("")
dev.off()


# Uniivariate models
model1 <- glm(ysnPatientDead ~ mdLMR, family=binomial(link='logit'),data=merge_sva_data)
summary(model)
exp(cbind(coef(model), confint(model)))

model2 <- glm(ysnPatientDead ~ mdNLR,family=binomial(link='logit'),data=merge_sva_data)
summary(model)
exp(cbind(coef(model), confint(model)))


# Multivariate models

model3 <- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
summary(model3)
exp(cbind(coef(model3), confint(model3)))

model4 <- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
summary(model4)
exp(cbind(coef(model4), confint(model4)))


#Lymphocyte to Monocyte Ratio (LMR) analysis
# 1. Age (Yes)
lm<-lm(mdLMR~Age+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data=merge_sva_data) #P= 0.000033 ***
summary(lm)

# 2. Gender (Yes)
lm<-lm(mdLMR~gender+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data=merge_sva_data)#P=0.00022 ***
summary(lm)


# 3. Drinking status (continous) (No)
lm<-lm(mdLMR~en_drink+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data=merge_sva_data)#P=0.45
summary(lm)

# 4. Smoking status (Ever/Never) (Yes)
lm<-lm(mdLMR~en_smoke+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data=merge_sva_data)#P=0.032 * 
summary(lm)

# 5. Tumour grade (Low/High) No
lm<-lm(mdLMR~grade+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data=merge_sva_data)#P=0.12 
summary(lm)

lm<-lm(mdLMR~grade+gender, data=data)
summary(lm)

lm<-lm(mdLMR~grade+gender+en_smoke+en_drink, data=data)
summary(lm)

lm<-lm(mdLMR~grade+gender+en_smoke+en_drink+Age, data=data)
summary(lm)

lm<-lm(mdLMR~as.factor(grade)+as.factor(gender)+as.factor(en_smoke)+en_drink+Age+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data=merge_sva_data)
summary(lm)

lm<-lm(mdLMR~as.factor(grade)+as.factor(gender)+as.factor(en_smoke)+en_drink+Age+as.factor(HPV16E6_pos2)+as.factor(SentrixID), data=data)
summary(lm)


# 6. HPV status (Yes)
lm<-lm(mdLMR~HPV16E6_pos2+as.factor(Sentrix_ID), data=data)
summary(lm)

cbind(coef(lm), confint(lm)) 

# 7. Barcode of samples (Did not work)

lm<-lm(mdLMR~IID, data=data)
summary(lm)

library(pROC)
# Comparing mdNLR and mdLMR for predicting death
fit1<- glm(ysnPatientDead ~ mdNLR+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
fit2<- glm(ysnPatientDead ~ mdLMR+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

tiff("ROC_HNSCC.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='green', print.auc = TRUE)
plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)
dev.off()

fit1<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
fit2<- glm(ysnPatientDead ~ mdNLR+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
fit3<- glm(ysnPatientDead ~ mdNLR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

tiff("ROC_HNSCC_multivariate.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='green', print.auc = TRUE)
plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = TRUE, print.auc.y = .2)
dev.off()




fit1<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
fit2<- glm(ysnPatientDead ~ mdLMR+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
fit3<- glm(ysnPatientDead ~ mdLMR+Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)

ci(roc1, method="bootstrap", parallel=TRUE)

preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)

tiff("ROC_HNSCC_multivariate_LMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='green', print.auc = TRUE)
plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)
plot(roc3, add=TRUE, col='blue', print.auc = TRUE, print.auc.y = .2)
dev.off()




fit1<- glm(ysnPatientDead ~ Age+as.factor(gender)+as.factor(grade)+as.factor(en_smoke)+as.factor(HPV16E6_pos2)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, family=binomial(link='logit'),data=merge_sva_data)
preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

tiff("ROC_HNSCC_multivariate.tif", res=300, compression = "lzw", height=5, width=10, units="in")
plot(roc1,col='green', print.auc = TRUE)
plot(roc2, add=TRUE, col='red', print.auc = TRUE, print.auc.y = .4)
dev.off()


#######################################
# myeloid cell subtypes in samples
######################################
myeloid<- read.csv("myeloid_final.csv", header=T)

box_myeloid <- myeloid[,c(24,2:6)]
names(box_myeloid)[1]<-paste("Overall_survival")

# plot the proportions by group of interest
colnames(box_myeloid)
box_myeloid$Overall_survival[box_myeloid$Overall_survival == 0] <- "Alive"
box_myeloid$Overall_survival[box_myeloid$Overall_survival == 1] <- "Dead"

stage.m <- melt(box_myeloid)
tiff("boxplot_myeloid.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Overall_survival)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Myeloid CpGs") +
  ylab("DNA methylation")+
  ggtitle("")
dev.off()


# Compare statistically if the difference is real
cg00901982_alv<-box_myeloid$cg00901982[box_myeloid$Overall_survival=="Alive"]
cg00901982_dead<-box_myeloid$cg00901982[box_myeloid$Overall_survival=="Dead"]
wilcox.test(cg00901982_alv, cg00901982_dead)# p-value = 0.0005
mean(cg00901982_alv)
sd(cg00901982_alv)
mean(cg00901982_dead)
sd(cg00901982_dead)


cg25938803_alv<-box_myeloid$cg25938803[box_myeloid$Overall_survival=="Alive"]
cg25938803_dead<-box_myeloid$cg25938803[box_myeloid$Overall_survival=="Dead"]
wilcox.test(cg25938803_alv, cg25938803_dead)# p-value = 0.0001
mean(cg25938803_alv)
sd(cg25938803_alv)
mean(cg25938803_dead)
sd(cg25938803_dead)


cg01591037_alv<-box_myeloid$cg01591037[box_myeloid$Overall_survival=="Alive"]
cg01591037_dead<-box_myeloid$cg01591037[box_myeloid$Overall_survival=="Dead"]
wilcox.test(cg01591037_alv, cg01591037_dead)# p-value = 0.003 
mean(cg01591037_alv)
sd(cg01591037_alv)
mean(cg01591037_dead)
sd(cg01591037_dead)


cg03621504_alv<-box_myeloid$cg03621504[box_myeloid$Overall_survival=="Alive"]
cg03621504_dead<-box_myeloid$cg03621504[box_myeloid$Overall_survival=="Dead"]
wilcox.test(cg03621504_alv, cg03621504_dead)# p-value = 0.01
mean(cg03621504_alv)
sd(cg03621504_alv)
mean(cg03621504_dead)
sd(cg03621504_dead)


cg10456459_alv<-box_myeloid$cg10456459[box_myeloid$Overall_survival=="Alive"]
cg10456459_dead<-box_myeloid$cg10456459[box_myeloid$Overall_survival=="Dead"]
wilcox.test(cg10456459_alv, cg10456459_dead)# p-value = 0.00007
mean(cg10456459_alv)
sd(cg10456459_alv)
mean(cg10456459_dead)
sd(cg10456459_dead)




## Disease free survival
box_myeloid <- myeloid[,c(30,2:6)]
names(box_myeloid)[1]<-paste("Disease_free_survival")

# plot the proportions by group of interest
colnames(box_myeloid)
box_myeloid$Disease_free_survival[box_myeloid$Disease_free_survival == 0] <- "Alive"
box_myeloid$Disease_free_survival[box_myeloid$Disease_free_survival == 1] <- "Dead"

# Remove NA values form the Disease Free Survival 
box_myeloid<-box_myeloid[complete.cases(box_myeloid), ]

stage.m <- melt(box_myeloid)
tiff("boxplot_myeloid_DFS.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Disease_free_survival)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 3) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Leucocyte subtype") +
  ylab("DNA methylation")+
  ggtitle("")
dev.off()


# Compare statistically if the difference is real
cg00901982_alv<-box_myeloid$cg00901982[box_myeloid$Disease_free_survival=="Alive"]
cg00901982_dead<-box_myeloid$cg00901982[box_myeloid$Disease_free_survival=="Dead"]
wilcox.test(cg00901982_alv, cg00901982_dead)# p-value = 0.0008




cg25938803_alv<-box_myeloid$cg25938803[box_myeloid$Disease_free_survival=="Alive"]
cg25938803_dead<-box_myeloid$cg25938803[box_myeloid$Disease_free_survival=="Dead"]
wilcox.test(cg25938803_alv, cg25938803_dead)# p-value = 0.02



cg01591037_alv<-box_myeloid$cg01591037[box_myeloid$Disease_free_survival=="Alive"]
cg01591037_dead<-box_myeloid$cg01591037[box_myeloid$Disease_free_survival=="Dead"]
wilcox.test(cg01591037_alv, cg01591037_dead)# p-value = 0.009


cg03621504_alv<-box_myeloid$cg03621504[box_myeloid$Disease_free_survival=="Alive"]
cg03621504_dead<-box_myeloid$cg03621504[box_myeloid$Disease_free_survival=="Dead"]
wilcox.test(cg03621504_alv, cg03621504_dead)# p-value = 0.03


cg10456459_alv<-box_myeloid$cg10456459[box_myeloid$Disease_free_survival=="Alive"]
cg10456459_dead<-box_myeloid$cg10456459[box_myeloid$Disease_free_survival=="Dead"]
wilcox.test(cg10456459_alv, cg10456459_dead)# p-value = 0.001


#######################################
# mdNLR mdLMR in HPV samples
######################################
box_HPV <- data[,c(17,10,11)]
names(box_HPV)[1]<-paste("HPV_Status")

# plot the proportions by group of interest
colnames(box_HPV)
box_HPV$HPV_Status[box_HPV$HPV_Status == 0] <- "Negative"
box_HPV$HPV_Status[box_HPV$HPV_Status == 1] <- "Positive"
library(ggplot2)
library(reshape2)
stage.m <- melt(box_HPV)
tiff("boxplot_HPV_mdNLR_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = HPV_Status)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 2) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("SI index measures") +
  ylab("SI index")+
  ggtitle("")
dev.off()

# Compare statistically if the difference is real (mdNLR)
mdNLR_hpvneg<-box_HPV$mdNLR[box_HPV$HPV_Status=="Negative"]
mdNLR_hpvpos<-box_HPV$mdNLR[box_HPV$HPV_Status=="Positive"]
wilcox.test(mdNLR_hpvneg, mdNLR_hpvpos)# p-value = 0.003

# Compare statistically if the difference is real (mdLMR)
mdLMR_hpvneg<-box_HPV$mdLMR[box_HPV$HPV_Status=="Negative"]
mdLMR_hpvpos<-box_HPV$mdLMR[box_HPV$HPV_Status=="Positive"]
wilcox.test(mdLMR_hpvneg, mdLMR_hpvpos)# p-value = 9.446e-05

#######################################
# mdNLR mdLMR in smoking samples
######################################
box_smoke <- data[,c(14,10,11)]
names(box_smoke)[1]<-paste("Smoking_Status")

# plot the proportions by group of interest
colnames(box_smoke)
box_smoke$Smoking_Status[box_smoke$Smoking_Status == 0] <- "Never"
box_smoke$Smoking_Status[box_smoke$Smoking_Status == 1] <- "Ever"
# many samples have NA to remove NA samples
box_smoke<-box_smoke[complete.cases(box_smoke), ] # 17 samples removed

library(ggplot2)
library(reshape2)
stage.m <- melt(box_smoke)
tiff("boxplot_smoke_mdNLR_mdLMR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(data = stage.m, aes(x = variable, y = value, fill = Smoking_Status)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 2) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("SI index measures") +
  ylab("SI index")+
  ggtitle("")
dev.off()

# Compare statistically if the difference is real (mdNLR)
mdNLR_never<-box_smoke$mdNLR[box_smoke$Smoking_Status=="Negative"]
mdNLR_ever<-box_smoke$mdNLR[box_smoke$Smoking_Status=="Positive"]
wilcox.test(mdNLR_never, mdNLR_ever)# p-value = 0.014

# Compare statistically if the difference is real (mdLMR)
mdLMR_never<-box_smoke$mdLMR[box_smoke$Smoking_Status=="Negative"]
mdLMR_ever<-box_smoke$mdLMR[box_smoke$Smoking_Status=="Positive"]
wilcox.test(mdLMR_never, mdLMR_ever)# p-value = 0.077





# Compare statistically if the difference is real
cg00901982_alv<-box_myeloid$cg00901982[box_myeloid$Overall_survival=="Alive"]
cg00901982_dead<-box_myeloid$cg00901982[box_myeloid$Overall_survival=="Dead"]
wilcox.test(cg00901982_alv, cg00901982_dead)# p-value = 0.0005

#####################################################################################

# new data analysis after obtaining file from Becky 07/03/18

#######################################################################################

setwd("//ads.bris.ac.uk/filestore/MyFiles/Staff16/sa16666/Documents/mDNLR project/HNSCC_mdNLR/H_N5000_V2.3")
library(readstata13)
dta<-read.dta13("H&N5000 V2.3.dta")

#############################################

# Merging the sva and dataframe 

#############################################
# Load data
load("~/mDNLR project/HNSCC_mdNLR/Data_Final_dataset_HNSCC_DFS.RData")
data<-new_df
rm(new_df)

library(tidyr)
unite<-unite(data,basename,Sentrix_ID,Sentrix_Pos, sep = "_")
sva<-read.table("sva_pheno_final.txt", header=T)
merge_sva_data<-merge(unite, sva, by.x = "basename", by.y = "me")

############################################################################################
# Merge merge_sva_data (basename) and dta (Epigenetics_ID) files
############################################################################################
total<-merge(merge_sva_data, dta, by.x = "basename", by.y = "Epigenetics_ID") # Only 415 samples were obatined
colnames(total)
write.csv(total, "total_merge_sva_HN5000data.csv")

###
# Remove columns which are not need
newdata<-total[,c(1:34,76,85,86)]

write.csv(newdata, "Data_analysis_HN5000data.csv")
# Change column names and call Alive as 1 and Dead as 2 in excel
#Load the dataset
newdata<-read.csv("Data_analysis_HN5000data.csv", header=T)

# Add values above and below median

median(newdata$mdNLR)#1.862476
median(newdata$mdLMR)#3.204831

newdata$mdNLR_ab_bel <- ifelse(newdata$mdNLR<= 1.862476, "Below Median",
                               ifelse(newdata$mdNLR> 1.862476, "Above Median",NA  ))

newdata$mdLMR_ab_bel <- ifelse(newdata$mdLMR<= 3.204831, "Below Median",
                               ifelse(newdata$mdLMR> 3.204831, "Above Median",NA  ))

###########################################
# Run cox models
############################################
library(survival)

# For mdNLR

# Univariate Model
res.cox <- coxph(Surv(Time, Survival_Status) ~ mdNLR_ab_bel, data = newdata) # P= 0.00264 **
summary(res.cox)
res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status)
res.cox <- coxph(Surv(Time, Survival_Status) ~ mdNLR_ab_bel+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(grade)+as.factor(HPV16E6_pos2.x), data = newdata) # P= 0.05917 .
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status and SVs)
res.cox <- coxph(Surv(Time, Survival_Status) ~ mdNLR_ab_bel+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(grade)+as.factor(HPV16E6_pos2.x)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data = newdata) # P= 0.07047 .
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# For mdLMR

# Univariate Model
res.cox <- coxph(Surv(Time, Survival_Status) ~ mdLMR_ab_bel, data = newdata) # P= 0.0000232 ***
summary(res.cox)
res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status)
res.cox <- coxph(Surv(Time, Survival_Status) ~ mdLMR_ab_bel+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(grade)+as.factor(HPV16E6_pos2.x), data = newdata) # P= 0.00139 ** 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status and SVs)
res.cox <- coxph(Surv(Time, Survival_Status) ~ mdLMR_ab_bel+Age+as.factor(gender)+as.factor(en_smoke)+as.factor(grade)+as.factor(HPV16E6_pos2.x)+sv.1+sv.2+sv.3+sv.4+sv.5+sv.6+sv.7+sv.8+sv.9+sv.10, data = newdata) # P= 0.000809 ***
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

################################################
# Kaplan Meier Curves
#################################################

library(survminer)
tiff("survplot_OS_mdNLR_HN5000.tif", res=300, compression = "lzw", height=8, width=15, units="in")
fit <- survfit(Surv(Time, Survival_Status) ~  mdNLR_ab_bel, data =  newdata)
ggsurvplot(fit, data = newdata, risk.table = TRUE)
dev.off()

tiff("survplot_OS_mdLMR_HN5000.tif", res=300, compression = "lzw", height=8, width=15, units="in")
fit <- survfit(Surv(Time, Survival_Status) ~ mdLMR_ab_bel, data = newdata)
ggsurvplot(fit, data = newdata, risk.table = TRUE)
dev.off()



setwd("//ads.bris.ac.uk/filestore/MyFiles/Staff16/sa16666/Documents/mDNLR project/HNSCC_mdNLR/New_ROC")
load("new_sva_1.Rdata")
dim(meth1)
table(summary(is.na(meth1)))

library(sva)


# The SVA model was not working so modified it by converting ysnPatientDead into 0= Alive and 1=Dead

setwd("//ads.bris.ac.uk/filestore/MyFiles/Staff16/sa16666/Documents/mDNLR project/HNSCC_mdNLR/New_ROC")

rm(pheno1)
pheno1<-read.csv("pheno1_HN5000_170318.csv", header=T)
colnames(pheno1)

mod<-model.matrix(~as.numeric(mdNLR)+as.numeric(mdLMR)+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(ysnPatient_C)+as.numeric(hn1_dv_death_cons)+as.factor(drinking), data=pheno1) 

mod0<-model.matrix(~as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(ysnPatient_C)+as.numeric(hn1_dv_death_cons)+as.factor(drinking), data=pheno1) 

sva<-sva(meth1,mod,mod0, n.sv=10)

sva<-as.data.frame(sva)

sva_1<-data.frame(sva$sv)

pheno_merge<-cbind(pheno1, sva_1)

names(pheno_merge)[82]<-paste("sv1")
names(pheno_merge)[83]<-paste("sv2")
names(pheno_merge)[84]<-paste("sv3")
names(pheno_merge)[85]<-paste("sv4")
names(pheno_merge)[86]<-paste("sv5")
names(pheno_merge)[87]<-paste("sv6")
names(pheno_merge)[88]<-paste("sv7")
names(pheno_merge)[89]<-paste("sv8")
names(pheno_merge)[90]<-paste("sv9")
names(pheno_merge)[91]<-paste("sv10")


# Add values above and below median

median(pheno_merge$mdNLR)#1.862476
median(pheno_merge$mdLMR)#3.204831

pheno_merge$mdNLR_ab_bel <- ifelse(pheno_merge$mdNLR<= 1.862476, "Below Median",
                               ifelse(pheno_merge$mdNLR> 1.862476, "Above Median",NA  ))

pheno_merge$mdLMR_ab_bel <- ifelse(pheno_merge$mdLMR<= 3.204831, "Below Median",
                               ifelse(pheno_merge$mdLMR> 3.204831, "Above Median",NA  ))

###########################################
# Run cox models
############################################
library(survival)
options(scipen=999)

# For mdNLR

# Univariate Model
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel, data = pheno_merge) # P=  0.00253 **
summary(res.cox)
res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model minimally adjusted (Age, Gender)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex), data = pheno_merge) # P=  0.012 * 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model minimally adjusted (Age, Gender, HPV)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(hpv), data = pheno_merge) # P=  0.02493 *  
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status, Drinking)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking), data = pheno_merge) # P= 0.046306 *
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status, Drinking and SVs)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data = pheno_merge) # P= 0.19706 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# Multivariate Model (Age, Gender, Grade, Smoking, HPV status, Drinking, ICD Site and SVs)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking)+as.factor(pheno1$hn1_ICD_group_conf)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data = pheno_merge) # P= 0.19706 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# For mdLMR

# Univariate Model
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel, data = pheno_merge) # P= 2.15e-05 ***
summary(res.cox)
res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model minimally adjusted (Age, Gender)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex), data = pheno_merge) # P=  0.000222 *** 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model minimally adjusted (Age, Gender, HPV)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(hpv), data = pheno_merge) # P=  0.000654 ***  
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# Multivariate Model (Age, Gender, Grade, Smoking, Drinking and HPV status)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking), data = pheno_merge) # P= 0.000848 ***
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, Drinking, HPV status and SVs)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data =pheno_merge) # P= 0.00188 ** 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# Multivariate Model (Age, Gender, Grade, Smoking, HPV status, Drinking, ICD Site and SVs)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking)+as.factor(pheno1$hn1_ICD_group_conf)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data = pheno_merge) # P= 0.19706 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


################################################
# Kaplan Meier Curves
#################################################

library(survminer)
tiff("survplot_OS_mdNLR_HN5000_170318.tif", res=300, compression = "lzw", height=8, width=15, units="in")
fit <- survfit(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel, data =  pheno_merge)
ggsurvplot(fit, data = pheno_merge, risk.table = TRUE)
dev.off()

tiff("survplot_OS_mdLMR_HN5000_170318.tif", res=300, compression = "lzw", height=8, width=15, units="in")
fit <- survfit(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel, data =  pheno_merge)
ggsurvplot(fit, data = pheno_merge, risk.table = TRUE)
dev.off()


###############################################
# Analysis of only OPC cases
###############################################
pheno_opc <- pheno1[which(pheno1$hn1_ICD_group_conf == "2 - Oropharynx"),]
dim(meth1)
colnames(meth1)
meth1_opc<-meth1[,colnames(meth1) %in% pheno_opc$me]


mod<-model.matrix(~as.numeric(mdNLR)+as.numeric(mdLMR)+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(ysnPatient_C)+as.numeric(hn1_dv_death_cons)+as.factor(drinking), data=pheno_opc) 

mod0<-model.matrix(~as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(ysnPatient_C)+as.numeric(hn1_dv_death_cons)+as.factor(drinking), data=pheno_opc) 

sva<-sva(meth1_opc,mod,mod0, n.sv=10)

sva_1<-data.frame(sva$sv)

pheno_merge<-cbind(pheno_opc, sva_1)

names(pheno_merge)[82]<-paste("sv1")
names(pheno_merge)[83]<-paste("sv2")
names(pheno_merge)[84]<-paste("sv3")
names(pheno_merge)[85]<-paste("sv4")
names(pheno_merge)[86]<-paste("sv5")
names(pheno_merge)[87]<-paste("sv6")
names(pheno_merge)[88]<-paste("sv7")
names(pheno_merge)[89]<-paste("sv8")
names(pheno_merge)[90]<-paste("sv9")
names(pheno_merge)[91]<-paste("sv10")


# Add values above and below median

median(pheno_merge$mdNLR)#1.836952
median(pheno_merge$mdLMR)#3.234055

pheno_merge$mdNLR_ab_bel <- ifelse(pheno_merge$mdNLR<= 1.836952, "Below Median",
                                   ifelse(pheno_merge$mdNLR> 3.234055, "Above Median",NA  ))

pheno_merge$mdLMR_ab_bel <- ifelse(pheno_merge$mdLMR<= 1.836952, "Below Median",
                                   ifelse(pheno_merge$mdLMR> 3.234055, "Above Median",NA  ))

###########################################
# Run cox models
############################################
library(survival)

# For mdNLR

# Univariate Model
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel, data = pheno_merge) # P=  0.0259 *
summary(res.cox)
res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model minimally adjusted (Age, Gender)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex), data = pheno_merge) # P=  0.0343 * 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model minimally adjusted (Age, Gender, HPV)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(hpv), data = pheno_merge) # P=  0.118   
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status, Drinking)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking), data = pheno_merge) # P= 0.195338 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status, Drinking and SVs)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdNLR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data = pheno_merge) # P= 0.290979 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# For mdLMR

# Univariate Model
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel, data = pheno_merge) # P= 0.00403 **
summary(res.cox)
res.cox
exp(cbind(HR=coef(res.cox), confint(res.cox)))

# Multivariate Model minimally adjusted (Age, Gender)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex), data = pheno_merge) # P=  0.00756 ** 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model minimally adjusted (Age, Gender, HPV)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(hpv), data = pheno_merge) # P=  0.0632 . 
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated


# Multivariate Model (Age, Gender, Grade, Smoking, HPV status and Drinking)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking), data = pheno_merge) # P=  0.128162
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

# Multivariate Model (Age, Gender, Grade, Smoking, HPV status and SVs)
res.cox <- coxph(Surv(hn1_dv_death_cons, ysnPatient_C) ~ mdLMR_ab_bel+as.numeric(age)+as.factor(sex)+as.factor(grade)+as.factor(en_smoke)+as.factor(hpv)+as.factor(drinking)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data =pheno_merge) # P= 0.034572 *
summary(res.cox)
# Test the Proportional Hazards Assumption of a Cox Regression
test.ph <- cox.zph(res.cox)
test.ph # Not violated

##############################
# Remake Supplementary Figure 1B
# Monocytes, Neutrophils and Lymphocytes
###############################
pheno_merge$Bcell <- as.numeric(as.character(pheno_merge$Bcell))
pheno_merge$CD4T <- as.numeric(as.character(pheno_merge$CD4T))
pheno_merge$CD8T <- as.numeric(as.character(pheno_merge$CD8T))
pheno_merge$Eos <- as.numeric(as.character(pheno_merge$Eos))
pheno_merge$NK <- as.numeric(as.character(pheno_merge$NK))

pheno_merge$lympho<-pheno_merge$Bcell+pheno_merge$CD4T+pheno_merge$CD8T+pheno_merge$Eos+pheno_merge$NK

test<-pheno_merge[,c(60,21,22,94)]
test$ysnPatient_C<-as.factor(test$ysnPatient_C)
head(test)

test$ysnPatient_C<-ifelse(test$ysnPatient_C==0,"Alive","Dead")

Mono_good<-test$Mono[test$ysnPatient_C=="Alive"]
Mono_poor<-test$Mono[test$ysnPatient_C=="Dead"]
wilcox.test(Mono_good,Mono_poor)#  0.00009888


Neu_good<-test$Neu[test$ysnPatient_C=="Alive"]
Neu_poor<-test$Neu[test$ysnPatient_C=="Dead"]
wilcox.test(Neu_good,Neu_poor)#  0.01951


lympho_good<-test$lympho[test$ysnPatient_C=="Alive"]
lympho_poor<-test$lympho[test$ysnPatient_C=="Dead"]
wilcox.test(lympho_good,lympho_poor)#  0.01951




library(reshape2)
test.m <- melt(test, id= "ysnPatient_C")
library(ggplot2)
tiff("boxplot_leukocytes_HNSCC5000.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = ysnPatient_C)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab(" ") +
  ylab(" ")+
  ggtitle("")
dev.off()

##############################################
# Figure 1C and 1D

#############################################
test<-pheno_merge[,c(60,24)]
test$ysnPatient_C<-as.factor(test$ysnPatient_C)
head(test)

test$ysnPatient_C<-ifelse(test$ysnPatient_C==0,"Alive","Dead")

res.by <- by(test$mdNLR,test$ysnPatient_C, mean)
res.by

mdNLR_good<-test$mdNLR[test$ysnPatient_C=="Alive"]
mdNLR_poor<-test$mdNLR[test$ysnPatient_C=="Dead"]
wilcox.test(mdNLR_good,mdNLR_poor)#  

library(reshape2)
test.m <- melt(test, id= "ysnPatient_C")
library(ggplot2)
tiff("boxplot_mdNLR_HNSCC5000.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = ysnPatient_C)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab(" ") +
  ylab(" ")+
  ggtitle("")
dev.off()


test<-pheno_merge[,c(60,25)]
test$ysnPatient_C<-as.factor(test$ysnPatient_C)
head(test)

test$ysnPatient_C<-ifelse(test$ysnPatient_C==0,"Alive","Dead")

res.by <- by(test$mdLMR,test$ysnPatient_C, mean)
res.by

mdLMR_good<-test$mdLMR[test$ysnPatient_C=="Alive"]
mdLMR_poor<-test$mdLMR[test$ysnPatient_C=="Dead"]
wilcox.test(mdLMR_good,mdLMR_poor)# 


library(reshape2)
test.m <- melt(test, id= "ysnPatient_C")
library(ggplot2)
tiff("boxplot_mdLMR_HNSCC5000.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = ysnPatient_C)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab(" ") +
  ylab(" ")+
  ggtitle("")
dev.off()


##############################################
# Supplementary Figure 3

#############################################
# HPV status
test<-pheno_merge[,c(15,24,25)]
test$HPV16E6_pos2.x<-as.factor(test$HPV16E6_pos2.x)
head(test)

test$HPV16E6_pos2.x<-ifelse(test$HPV16E6_pos2.x==0,"Negative","Positive")

library(reshape2)
test.m <- melt(test, id= "HPV16E6_pos2.x")
library(ggplot2)
tiff("boxplot_mdNLR_mdLMR_HPV_HNSCC5000.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = HPV16E6_pos2.x)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab(" ") +
  ylab(" ")+
  ggtitle("")
dev.off()

# Smoking status

test<-pheno_merge[,c(9,24,25)]
test$evernever<-as.factor(test$evernever)
head(test)

test$evernever<-ifelse(test$evernever==0,"Never","Ever")

library(reshape2)
test.m <- melt(test, id= "evernever")
library(ggplot2)
tiff("boxplot_mdNLR_mdLMR_SMK_HNSCC5000.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = evernever)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab(" ") +
  ylab(" ")+
  ggtitle("")
dev.off()

###################################################
# Remkaing Table4 and Supplementary figure 4 
###################################################
# Load myeloid final recieved from Ryan earlier
myeloid<-read.csv("myeloid_final.csv", header=T)
# Select only columns of interest
colnames(myeloid)
myeloid_sel<-myeloid[,c(1:6)]

# Merge the pheno_merge data with myeloid sel data

data <- merge(pheno_merge, myeloid_sel, by.x = "me", by.y = "Barcode", all.x = T)
colnames(data)

data_dead<-data[data$ysnPatient_C == "1",]
data_alive<-data[data$ysnPatient_C == "0",]


res.by <- by(data$cg25938803, data$ysnPatient_C, mean)
res.by
res.by <- by(data$cg25938803, data$ysnPatient_C, sd)
res.by
options(scipen=999)
cg25938803_good<-data$cg25938803[data$ysnPatient_C=="0"]
cg25938803_poor<-data$cg25938803[data$ysnPatient_C=="1"]
wilcox.test(cg25938803_good,cg25938803_poor)#  0.001154


res.by <- by(data$cg10456459, data$ysnPatient_C, mean)
res.by
res.by <- by(data$cg10456459, data$ysnPatient_C, sd)
res.by
options(scipen=999)
cg10456459_good<-data$cg10456459[data$ysnPatient_C=="0"]
cg10456459_poor<-data$cg10456459[data$ysnPatient_C=="1"]
wilcox.test(cg10456459_good,cg10456459_poor)#  0.0003076


res.by <- by(data$cg01591037, data$ysnPatient_C, mean)
res.by
res.by <- by(data$cg01591037, data$ysnPatient_C, sd)
res.by
options(scipen=999)
cg01591037_good<-data$cg01591037[data$ysnPatient_C=="0"]
cg01591037_poor<-data$cg01591037[data$ysnPatient_C=="1"]
wilcox.test(cg01591037_good,cg01591037_poor)#  0.02488


res.by <- by(data$cg03621504, data$ysnPatient_C, mean)
res.by
res.by <- by(data$cg03621504, data$ysnPatient_C, sd)
res.by
options(scipen=999)
cg03621504_good<-data$cg03621504[data$ysnPatient_C=="0"]
cg03621504_poor<-data$cg03621504[data$ysnPatient_C=="1"]
wilcox.test(cg03621504_good,cg03621504_poor)#  0.03133


res.by <- by(data$cg00901982, data$ysnPatient_C, mean)
res.by
res.by <- by(data$cg00901982, data$ysnPatient_C, sd)
res.by
options(scipen=999)
cg00901982_good<-data$cg00901982[data$ysnPatient_C=="0"]
cg00901982_poor<-data$cg00901982[data$ysnPatient_C=="1"]
wilcox.test(cg00901982_good,cg00901982_poor)#  0.0006635


#################################################
# Boxplot
# myeloid CpGs
#################################################

test<-data[,c(60,95:99)]
test$ysnPatient_C<-as.factor(test$ysnPatient_C)
head(test)

test$ysnPatient_C<-ifelse(test$ysnPatient_C==0,"Alive","Dead")

library(reshape2)
test.m <- melt(test, id= "ysnPatient_C")
library(ggplot2)
tiff("boxplot_myeloid_CPG_HNSCC5000.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = ysnPatient_C)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  #scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab(" ") +
  ylab(" ")+
  ggtitle("")
dev.off()




