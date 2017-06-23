# Epigenetic Age for Cleft study - 14/11/2016 - Gemma Sharp #

# Raw, non-normalised beta values for all samples uploaded to https://dnamage.genetics.ucla.edu plus an annotation file containing chronological age, sex and tissue
# Normalisation and advanced analysis for blood options selected
# Output file: Full_Aries_for_Horvath.output 2.csv

load("Cleft/NoControls/betas.raw.Robj")
load("Cleft/NoControls/samplesheet.cleft.Robj")
samplesheet.cleft<-samplesheet.cleft[which(samplesheet.cleft$sample_type=="Blood"),]
samplesInBoth<-intersect(colnames(betas.raw),samplesheet.cleft$Sample_Name)
betas.raw<-betas.raw[,samplesInBoth]
samplesheet.cleft<-samplesheet.cleft[match(samplesInBoth,samplesheet.cleft$Sample_Name),]
all(colnames(betas.raw)==samplesheet.cleft$Sample_Name)

dat0=data.frame(betas.raw)
dat0=cbind(rownames(dat0), dat0)
datMiniAnnotation=read.csv("Cleft/datMiniAnnotation.csv")
match1=match(datMiniAnnotation[,1], dat0[,1] )
dat0Reduced=dat0[match1,]
dat0Reduced[,1]=as.character(dat0Reduced[,1])
dat0Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
datout=data.frame(dat0Reduced)
# make sure you output numeric variables...
for (i in 2:dim(datout)[[2]] ){datout[,i]= as.numeric(as.character(gsub(x=datout[,i],pattern="\"",replacement=""))) }
#replace "MethylationData" with a filename of your choice
write.table(datout,"Cleft/NoControls/MethylationDataForHorvath.csv", row.names=F, sep="," )
samplesheet.cleft$Age<-samplesheet.cleft$c_age/52 #age in years
samplesheet.cleft$Tissue<-"Blood"
samplesheet.cleft$Female[samplesheet.cleft$Sex=="F"]<-1
samplesheet.cleft$Female[samplesheet.cleft$Sex=="M"]<-0
write.table(samplesheet.cleft,"Cleft/NoControls/PhenoDataForHorvath.csv", row.names=F, sep=",")

agedat<-read.csv("/panfs/panasas01/sscm/gs8094/Cleft/NoControls/MethylationDataForHorvath.output.csv",stringsAsFactors=F)

#Horvath suggests (in https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/TUTORIALonlineCalculator.pdf):
#To assess whether DNAmAge relates to a disease outcome, I use the following covariate list DNAmAge+Age+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran. Obviously, you would also adjust for standard variables such as gender, race, body mass index, prior history of disease e.g. cancer, type II diabetes status, etc. 

cor(agedat$DNAmAge,agedat$Age,use="pairwise")#0.7 (spearman 0.6)
cor(agedat[which(agedat$cleft_type=="CPO"),"DNAmAge"],agedat[which(agedat$cleft_type=="CPO"),"Age"],use="pairwise")#0.6 (spearman 0.5)
cor(agedat[which(agedat$cleft_type=="CLO"),"DNAmAge"],agedat[which(agedat$cleft_type=="CLO"),"Age"],use="pairwise")#0.3 (spearman 0.3)
cor(agedat[which(agedat$cleft_type=="CLP"),"DNAmAge"],agedat[which(agedat$cleft_type=="CLP"),"Age"],use="pairwise")#0.7 (spearman 0.3)


summary(lm(CLO.CLP~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#                   Estimate Std. Error t value Pr(>|t|)
#DNAmAge          -0.2069970  0.1552108  -1.334   0.1859 
#Age              -0.6447581  0.3858741  -1.671   0.0985 
summary(lm(CPO.CLP~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#                   Estimate Std. Error t value Pr(>|t|)
#DNAmAge          -0.1179281  0.1037743  -1.136    0.259    
#Age               1.2106810  0.2072050   5.843 1.06e-07
summary(lm(CPO.CLO~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#                   Estimate Std. Error t value Pr(>|t|)
#DNAmAge          -0.0370836  0.0765859  -0.484    0.630    
#Age               1.3776139  0.1373004  10.034  9.5e-16


summary(lm(smoking.score.reese.blood~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#                   Estimate Std. Error t value Pr(>|t|)
#DNAmAge          -4.424e-02  1.001e-01  -0.442   0.6591  
#Age              -1.675e-01  1.887e-01  -0.888   0.3764 
summary(lm(m_age~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#                   Estimate Std. Error t value Pr(>|t|)
#DNAmAge          -0.275711   1.249953  -0.221 0.825776    
#Age               1.080759   2.357757   0.458 0.647460
summary(lm(parity>1~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#                   Estimate Std. Error t value Pr(>|t|)
#DNAmAge           0.0110288  0.1655930   0.067  0.94710   
#Age              -0.4575417  0.3372955  -1.357  0.17949
summary(lm(m_smoking_conception~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))
#DNAmAge           0.157707   0.198193   0.796    0.432
#Age               0.155058   0.436321   0.355    0.725 
summary(lm(m_smoking_conception~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat))

agedat$matedu2<-NA
agedat$matocc2<-NA
agedat$matedu2[which(agedat$m_education %in% c(9,10))]<-"Uni"
agedat$matedu2[which(agedat$m_education %in% c(1,2,3,4,5,6,7,8,11,12,13,15))]<-"NoUni"
agedat$matocc2[which(agedat$m_occupation %in% c(1,2,3))]<-"non-manual"
agedat$matocc2[which(agedat$m_occupation %in% c(4,5,6,7,8,9,10))]<-"manual or unskilled"

summary(glm(as.factor(matedu2)~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat,family="binomial"))
#DNAmAge          -0.600925   0.731068  -0.822    0.411
#Age               1.916519   1.497937   1.279    0.201

summary(glm(as.factor(matocc2)~DNAmAge+Age+Sex+CD8.naive + CD8pCD28nCD45RAn + PlasmaBlast+CD4T+NK+Mono+Gran,data=agedat,family="binomial"))
#DNAmAge           0.320525   0.760384   0.422    0.673
#Age               1.721325   1.574430   1.093    0.274

summary(agedat$DNAmAge*52)
#Age was predicted to be between 20 weeks before birth to 3 years old, with a median of 36 weeks (IQR 23 to 52)
#Actual age ranged from 10 weeks to 85 weeks (21 months. 1.6 years), with a medium of 20 weeks (IQR 15 to 38)
#Deviations between actual and predicted age were not associated with cleft subtype or any measured covariates (after adjustment for estimated cell counts and sex)

require(data.table)
agedat<-data.table(agedat)
agedat$DNAmAge_mnths<-agedat$DNAmAge*12
agedat[,list(Mean=mean(DNAmAge_mnths,na.rm=T),Sd=sd(DNAmAge_mnths,na.rm=T),LCI=mean(DNAmAge_mnths,na.rm=T)-1.96*(sd(DNAmAge_mnths,na.rm=T)/sqrt(length(DNAmAge_mnths))),UCI=mean(DNAmAge_mnths,na.rm=T)+1.96*(sd(DNAmAge_mnths,na.rm=T)/sqrt(length(DNAmAge_mnths)))),by=cleft_type]
summary(aov(DNAmAge_mnths~cleft_type,data=agedat))

agedat$AA<-agedat$AAHOAdjCellCounts*12
agedat[,list(Mean=mean(AA,na.rm=T),Sd=sd(AA,na.rm=T),LCI=mean(AA,na.rm=T)-1.96*(sd(AA,na.rm=T)/sqrt(length(AA))),UCI=mean(AA,na.rm=T)+1.96*(sd(AA,na.rm=T)/sqrt(length(AA)))),by=cleft_type]
summary(aov(AA~cleft_type,data=agedat))

t.test(AA~CLO.CLP,data=agedat)$p.value #0.13
t.test(AA~CPO.CLO,data=agedat)$p.value #0.85
t.test(AA~CPO.CLP,data=agedat)$p.value #0.16




