#Inferring maternal smoking status from cleft blood methylation data

setwd("Cleft/NoControls/")
load("norm.beta.pc15.Robj")
load("samplesheet.cleft.Robj")
samplesheet<-samplesheet.cleft
bloodsamples<-samplesheet[match(colnames(norm.beta),samplesheet$Sample_Name),]
bloodsamples<-bloodsamples[which(bloodsamples$sample_type=="Blood"),]
methylation.blood<-norm.beta[,match(bloodsamples$Sample_Name,colnames(norm.beta))]
tissuesamples<-samplesheet[match(colnames(norm.beta),samplesheet$Sample_Name),]
tissuesamples<-tissuesamples[which(tissuesamples$sample_type%in%c("Lip","Palate")),]
methylation.tissue<-norm.beta[,match(tissuesamples$Sample_Name,colnames(norm.beta))]

pace <- read.csv("/panfs/panasas01/sscm/gs8094/Cleft/pace-prenatal-smoking-associations.csv",stringsAsFactors=F)
reese <- read.csv("/panfs/panasas01/sscm/gs8094/Cleft/Reese.csv",stringsAsFactors=F)

library(mixtools)
library(matrixStats)
library(limma)

colnames(pace) <- tolower(colnames(pace))
pace<-na.omit(pace[match(row.names(methylation.blood),pace$cpg),])
reese<-na.omit(reese[match(row.names(methylation.blood),reese$cpg),])

meth.score <- function(methylation, weights, reference) {
    cpg.idx <- match(reference$cpg, rownames(methylation))
    scale(t(methylation[cpg.idx,,drop=F])) %*% weights
}

pace.smoking.coefs <- pace$coef
names(pace.smoking.coefs) <- pace$cpg
pace.cpg.weights <- pace.smoking.coefs/sum(abs(pace.smoking.coefs))
prediction.pace.blood <- meth.score(methylation.blood, pace.cpg.weights,pace)
prediction.pace.tissue <- meth.score(methylation.tissue, pace.cpg.weights,pace)

reese.smoking.coefs <- reese$coef
names(reese.smoking.coefs) <- reese$cpg
reese.cpg.weights <- reese.smoking.coefs/sum(abs(reese.smoking.coefs))
prediction.reese.blood <- meth.score(methylation.blood, reese.cpg.weights,reese)
prediction.reese.tissue <- meth.score(methylation.tissue, reese.cpg.weights,reese)

remove.pcs <- function(methylation) {
    var.idx <- order(rowVars(methylation, na.rm=T), decreasing=T)
    fit <- prcomp(t(methylation[var.idx[1:10000],]), scale=T)
    fit <- lmFit(methylation, design=fit$x[,1:10])
    methylation.adjusted <- residuals(fit, methylation)
}

methylation.blood.adj <- remove.pcs(methylation.blood)
methylation.tissue.adj <- remove.pcs(methylation.tissue)

prediction.pace.blood.adj <- meth.score(methylation.blood.adj, pace.cpg.weights,pace)
prediction.pace.tissue.adj <- meth.score(methylation.tissue.adj, pace.cpg.weights,pace)
prediction.reese.blood.adj <- meth.score(methylation.blood.adj, reese.cpg.weights,reese)
prediction.reese.tissue.adj <- meth.score(methylation.tissue.adj, reese.cpg.weights,reese)

#fit <- normalmixEM(prediction)
#prediction.bin <- apply(fit$posterior, 1, which.max) -1
#if (fit$mu[1] > fit$mu[2])
#    prediction.bin <- 1-prediction.bin

#fit.adj <- normalmixEM(prediction.adj)
#prediction.adj.bin <- apply(fit.adj$posterior, 1, which.max)-1
#if (fit.adj$mu[1] > fit.adj$mu[2])
#    prediction.adj.bin <- 1-prediction.adj.bin

Blood<-merge(samplesheet[,c("Sample_Name","participant","sample_type")],prediction.pace.blood,by.x="Sample_Name",by.y="row.names",all=T)
Blood<-merge(Blood,prediction.reese.blood,by.x="Sample_Name",by.y="row.names",all=T)
Blood<-merge(Blood,prediction.pace.blood.adj,by.x="Sample_Name",by.y="row.names",all=T)
Blood<-merge(Blood,prediction.reese.blood.adj,by.x="Sample_Name",by.y="row.names",all=T)
colnames(Blood)<-c("Sample_Name","participant","sample_type","smoking.score.pace.blood",
"smoking.score.reese.blood","smoking.score.pace.blood.adj","smoking.score.reese.blood.adj")

Tissue<-merge(samplesheet[,c("Sample_Name","participant","sample_type")],prediction.pace.tissue,by.x="Sample_Name",by.y="row.names",all=T)
Tissue<-merge(Tissue,prediction.reese.tissue,by.x="Sample_Name",by.y="row.names",all=T)
Tissue<-merge(Tissue,prediction.pace.tissue.adj,by.x="Sample_Name",by.y="row.names",all=T)
Tissue<-merge(Tissue,prediction.reese.tissue.adj,by.x="Sample_Name",by.y="row.names",all=T)
colnames(Tissue)<-c("Sample_Name","participant","sample_type","smoking.score.pace.tissue",
"smoking.score.reese.tissue","smoking.score.pace.tissue.adj","smoking.score.reese.tissue.adj")

Unique<-unique(as.character(samplesheet$participant))
Blood<-na.omit(Blood)
Tissue<-na.omit(Tissue)
All<-merge(Blood[,-c(1,3)],Tissue[,-c(1,3)],by="participant",all=T)
All$participant<-as.character(All$participant)
All->All.norm
save(All,file="Smoking_Status_Score.Rdata")







#######################################################################################
#### Testing score options


### RAW METHYLATION ###

#Inferring maternal smoking status from cleft blood methylation data

setwd("Cleft/NoControls/")
load("betas.raw.Robj")
load("samplesheet.cleft.Robj")
samplesheet<-samplesheet.cleft
bloodsamples<-samplesheet[match(colnames(betas.raw),samplesheet$Sample_Name),]
bloodsamples<-bloodsamples[which(bloodsamples$sample_type=="Blood"),]
methylation.blood<-betas.raw[,match(bloodsamples$Sample_Name,colnames(betas.raw))]
tissuesamples<-samplesheet[match(colnames(betas.raw),samplesheet$Sample_Name),]
tissuesamples<-tissuesamples[which(tissuesamples$sample_type%in%c("Lip","Palate")),]
methylation.tissue<-betas.raw[,match(tissuesamples$Sample_Name,colnames(betas.raw))]

pace <- read.csv("/panfs/panasas01/sscm/gs8094/Cleft/pace-prenatal-smoking-associations.csv",stringsAsFactors=F)
reese <- read.csv("/panfs/panasas01/sscm/gs8094/Cleft/Reese.csv",stringsAsFactors=F)

library(mixtools)
library(matrixStats)
library(limma)

colnames(pace) <- tolower(colnames(pace))
pace<-na.omit(pace[match(row.names(methylation.blood),pace$cpg),])
reese<-na.omit(reese[match(row.names(methylation.blood),reese$cpg),])

meth.score <- function(methylation, weights, reference) {
    cpg.idx <- match(reference$cpg, rownames(methylation))
    scale(t(methylation[cpg.idx,,drop=F])) %*% weights
}

pace.smoking.coefs <- pace$coef
names(pace.smoking.coefs) <- pace$cpg
pace.cpg.weights <- pace.smoking.coefs/sum(abs(pace.smoking.coefs))
prediction.pace.blood <- meth.score(methylation.blood, pace.cpg.weights,pace)
prediction.pace.tissue <- meth.score(methylation.tissue, pace.cpg.weights,pace)

reese.smoking.coefs <- reese$coef
names(reese.smoking.coefs) <- reese$cpg
reese.cpg.weights <- reese.smoking.coefs/sum(abs(reese.smoking.coefs))
prediction.reese.blood <- meth.score(methylation.blood, reese.cpg.weights,reese)
prediction.reese.tissue <- meth.score(methylation.tissue, reese.cpg.weights,reese)

remove.pcs <- function(methylation) {
    var.idx <- order(rowVars(methylation, na.rm=T), decreasing=T)
    fit <- prcomp(t(methylation[var.idx[1:10000],]), scale=T)
    fit <- lmFit(methylation, design=fit$x[,1:10])
    methylation.adjusted <- residuals(fit, methylation)
}

methylation.blood.adj <- remove.pcs(methylation.blood)
methylation.tissue.adj <- remove.pcs(methylation.tissue)

prediction.pace.blood.adj <- meth.score(methylation.blood.adj, pace.cpg.weights,pace)
prediction.pace.tissue.adj <- meth.score(methylation.tissue.adj, pace.cpg.weights,pace)
prediction.reese.blood.adj <- meth.score(methylation.blood.adj, reese.cpg.weights,reese)
prediction.reese.tissue.adj <- meth.score(methylation.tissue.adj, reese.cpg.weights,reese)

#fit <- normalmixEM(prediction)
#prediction.bin <- apply(fit$posterior, 1, which.max) -1
#if (fit$mu[1] > fit$mu[2])
#    prediction.bin <- 1-prediction.bin

#fit.adj <- normalmixEM(prediction.adj)
#prediction.adj.bin <- apply(fit.adj$posterior, 1, which.max)-1
#if (fit.adj$mu[1] > fit.adj$mu[2])
#    prediction.adj.bin <- 1-prediction.adj.bin

Blood<-merge(samplesheet[,c("Sample_Name","participant","sample_type")],prediction.pace.blood,by.x="Sample_Name",by.y="row.names",all=T)
Blood<-merge(Blood,prediction.reese.blood,by.x="Sample_Name",by.y="row.names",all=T)
Blood<-merge(Blood,prediction.pace.blood.adj,by.x="Sample_Name",by.y="row.names",all=T)
Blood<-merge(Blood,prediction.reese.blood.adj,by.x="Sample_Name",by.y="row.names",all=T)
colnames(Blood)<-c("Sample_Name","participant","sample_type","smoking.score.pace.blood",
"smoking.score.reese.blood","smoking.score.pace.blood.adj","smoking.score.reese.blood.adj")

Tissue<-merge(samplesheet[,c("Sample_Name","participant","sample_type")],prediction.pace.tissue,by.x="Sample_Name",by.y="row.names",all=T)
Tissue<-merge(Tissue,prediction.reese.tissue,by.x="Sample_Name",by.y="row.names",all=T)
Tissue<-merge(Tissue,prediction.pace.tissue.adj,by.x="Sample_Name",by.y="row.names",all=T)
Tissue<-merge(Tissue,prediction.reese.tissue.adj,by.x="Sample_Name",by.y="row.names",all=T)
colnames(Tissue)<-c("Sample_Name","participant","sample_type","smoking.score.pace.tissue",
"smoking.score.reese.tissue","smoking.score.pace.tissue.adj","smoking.score.reese.tissue.adj")

Unique<-unique(as.character(samplesheet$participant))
Blood<-na.omit(Blood)
Tissue<-na.omit(Tissue)
All<-merge(Blood[,-c(1,3)],Tissue[,-c(1,3)],by="participant",all=T)
All$participant<-as.character(All$participant)
All->All.raw
save(All,file="Smoking_Status_Score.Rdata")

#
The correlation between the score generated using normalised and raw methylation values is >0.9 for all methods, so I have no preference over which to use. The normalisation technique isn't doing anything to the data that makes it difficult to predict smoking.

There is weak correlation between the scores generated using tissue and blood.

cor(All.norm$smoking.score.reese.blood,All.norm$smoking.score.reese.tissue,use="pairwise")
#0.1692306
cor(All.norm$smoking.score.pace.blood,All.norm$smoking.score.pace.tissue,use="pairwise")
#0.1858598
cor(All.norm$smoking.score.reese.blood.adj,All.norm$smoking.score.reese.tissue.adj,use="pairwise")
#0.1562681
cor(All.norm$smoking.score.pace.blood.adj,All.norm$smoking.score.pace.tissue.adj,use="pairwise")
#0.1782736

#Smoking at conception isn't associated with any of the scores. PACE does better than REESE, but adjusting PACE for PCs reverses sign of association.
#The only exception is reese for tissue, which is associated with smoking at conception but in the wrong direction (smokers have lower score).

summary(lm(samplesheet$smoking.score.pace.blood~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)
(Intercept)                      0.008135   0.025878   0.314    0.754
samplesheet$m_smoking_conception 0.055728   0.051165   1.089    0.279
summary(lm(samplesheet$smoking.score.pace.blood.adj~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)
(Intercept)                       0.01440    0.01274   1.130    0.262
samplesheet$m_smoking_conception -0.02549    0.02518  -1.012    0.314
summary(lm(samplesheet$smoking.score.reese.blood~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)  
(Intercept)                      -0.12155    0.04837  -2.513   0.0139 *
samplesheet$m_smoking_conception  0.01706    0.09564   0.178   0.8589 
summary(lm(samplesheet$smoking.score.reese.blood.adj~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)  
(Intercept)                      -0.10523    0.04564  -2.306   0.0236 *
samplesheet$m_smoking_conception  0.01344    0.09023   0.149   0.8820  

summary(lm(samplesheet$smoking.score.pace.tissue~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)
(Intercept)                      -0.00993    0.01059  -0.938    0.351
samplesheet$m_smoking_conception  0.02724    0.02142   1.272    0.207
summary(lm(samplesheet$smoking.score.pace.tissue.adj~samplesheet$m_smoking_conception))
                                  Estimate Std. Error t value Pr(>|t|)
(Intercept)                      -0.008727   0.005994  -1.456    0.149
samplesheet$m_smoking_conception  0.001772   0.012124   0.146    0.884
summary(lm(samplesheet$smoking.score.reese.tissue~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                       0.13261    0.05047   2.627 0.010152 *  
samplesheet$m_smoking_conception -0.40837    0.10209  -4.000 0.000132 ***
summary(lm(samplesheet$smoking.score.reese.tissue.adj~samplesheet$m_smoking_conception))
                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                       0.12174    0.04921   2.474 0.015276 *  
samplesheet$m_smoking_conception -0.39769    0.09952  -3.996 0.000134 ***

#Reese score in blood is associated with missingness
samplesheet$smoking_missing<-0
samplesheet$smoking_missing[is.na(samplesheet$m_smoking_conception)]<-1

summary(lm(samplesheet$smoking.score.reese.blood~samplesheet$smoking_missing)) 
                            Estimate Std. Error t value Pr(>|t|)   
(Intercept)                 -0.11719    0.04370  -2.682  0.00774 **
samplesheet$smoking_missing  0.16564    0.05195   3.188  0.00159 **

#Blood EWAS of smoking at conception:
load("Cleft/NoControls/norm.beta.pc15.IQRtrimmed.Robj")
load("Cleft/NoControls/samplesheet.cleft.Robj")
samplesheet.cleft->samplesheet
norm.beta->meth

samplesheet.blood <- droplevels(samplesheet[which(samplesheet$sample_type=="Blood"),])
samplesheet.lip <- droplevels(samplesheet[which(samplesheet$sample_type=="Lip"),])
samplesheet.palate <- droplevels(samplesheet[which(samplesheet$sample_type=="Palate"),])
samplesheet.tissue <- droplevels(rbind(samplesheet.lip,samplesheet.palate))

require(CpGassoc)

load("/panfs/panasas01/sscm/gs8094/Common_files/fdata_new.RData")
noXY <- as.character(fdata.new[-which(fdata.new$CHR %in% c("X", "Y")), "TargetID"])
noXYSNP <- as.character(noXY[grep("rs", noXY, invert=TRUE)])
filtered<- fdata.new[fdata.new$TargetID %in% noXYSNP,]
meth <- subset(meth, row.names(meth) %in% filtered$TargetID)

dat1<-samplesheet.blood

#match betas and pheno data
meth1 <- meth[, colnames(meth) %in% dat1$Sample_Name]
dat1 <- dat1[match(colnames(meth1), dat1$Sample_Name), , drop=FALSE]
stopifnot(all(dat1$Sample_Name == colnames(meth1)))

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$m_smoking_conception))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->m_smoking_conception #Top hits not smoking related

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.reese.blood))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->reese_blood # Top hits are smoking related

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.reese.blood.adj))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->reese_blood_adj # Top hits are smoking related

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.pace.blood))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->pace_blood #Top hits not smoking related and lots of inflation

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.pace.blood.adj))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->pace_blood_adj #Top hits not smoking related and lots of inflation

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.reese.tissue))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->reese_tissue #Top hits not smoking related

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.reese.tissue.adj))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->reese_tissue_adj #Top hits not smoking related

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.pace.tissue))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->pace_tissue

res<-cpg.assoc(beta.val=meth1,indep=as.numeric(dat1$smoking.score.pace.tissue.adj))
res<-merge(res$results,res$coefficients,by.x="CPG.Labels",by.y="row.names",all=T)
res<-merge(res,fdata.new,by.x="CPG.Labels",by.y="TargetID",all.x=T,all.y=F)
res<-res[order(res$P.value),]
res->pace_tissue_adj

#coefficients from Reese correlate with coefficients from reese_blood EWAS
RB<-merge(reese_blood[,c("CPG.Labels","effect.size","P.value")],reese,by.x="CPG.Labels",by.y="cpg",all.x=F,all.y=T)
RBA<-merge(reese_blood_adj[,c("CPG.Labels","effect.size","P.value")],reese,by.x="CPG.Labels",by.y="cpg",all.x=F,all.y=T)

> cor(RB$effect.size,RB$coefficient)
[1] 0.5255082
> cor(RBA$effect.size,RBA$coefficient)
[1] 0.5570357

#Slightly better correlation when reese_blood is adjusted for PCs, so use reese_blood_adj

