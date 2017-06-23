# Code to QC and normalise the Cleft Collective 450k data using meffil #
# NOW WITH ADDED CONTROLS! #
# May 8th 2017 #
# Gemma Sharp #

# Would not run interactively using qsub -I -l walltime=12:00:00 as could not access RDSF data via qsub #
# Have to run on login node until find better alternative #

#######################################################################################################################################

######################
# Create samplesheet #
######################

setwd("Cleft/Controls2")
library(meffil)

#Basenames (Slides is the first group that were released on Oct 12, Slides2 is the second group released on Oct 21)
Slides<- c("200687170051", "200687170098", "200687170135", "200693590007","200693590015","200693590022","200693590023","200693590024", "200693590025","200693590028", "200693590032","200693590033","200693590037","200693590038","200693590060","200693590104")
Slides2<- c("200687170081",  "200687170158",  "200687170160",  "200693590035",  "200693590105","200687170143",  "200687170159", "200693590026",  "200693590056")
Basenames<-paste0("/projects/Cleft_Collective/Cleft/",Slides)
Basenames2<-paste0("/projects/Cleft_Collective/Cleft-Extra-2016-10-21/",Slides2)
Basenames<-c(Basenames,Basenames2)
samplesheet.cleft <- do.call("rbind",lapply(Basenames,meffil.create.samplesheet))

#Combine samplesheet with other data we have on these samples
load("/panfs/panasas01/sscm/gs8094/Cleft/phen.cleaned.Robj") #this was created from the raw pheno data using the script cleaning_pheno_450k_CC.R
samplesheet.cleft$Sex<-NULL
samplesheet.cleft<-merge(samplesheet.cleft,phen,by.x="Sample_Name",by.y="sentrix_id",all.x=T)
samplesheet.cleft$Study<-"CCC"
samplesheet.cleft$sample_type<-as.character(samplesheet.cleft$sample_type)
save(samplesheet.cleft,file="samplesheet.cleft.Robj")

############################################################
# Downloading and creating a samplesheet for some controls #
############################################################

# Possible controls for cleft children blood samples age 10 to 85 weeks, median 20

#GSE72556 age 3 to 5 Latino saliva (confounded by batch, tissue and ethnicity)
#GSE36054 GSE36064 (same kids) age 12 months to 203 months PBL (confounded by batch and tissue)
#GSE67444 Dried blood spots 6 to 60 months (confounded by batch, tissue, not sure on ethnicity)
#GSE62219 PBLs age 3 to 60 months all girls (confounded by batch, tissue, not sure about ethnicity)
#GSE59592 3 to 6 months blood Gambia (confounded by ethnicity, batch) raw data not available
#GSE54399 venous whole blood not sure of age (confounded by batch) raw data not available

# Downloading idats from GEO (if available)
# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery")
# require(GEOquery)

#GSE54399 <- getGEOSuppFiles("GSE54399") #None available
#untar("panfs/panasas01/sscm/gs8094/Cleft/GSE54399/GSE54399_RAW.tar", exdir = "panfs/panasas01/sscm/gs8094/Cleft/GSE54399/idat")

#GSE59592<- getGEOSuppFiles("GSE59592") #None available
#untar("panfs/panasas01/sscm/gs8094/Cleft/GSE59592/GSE59592_RAW.tar", exdir = "panfs/panasas01/sscm/gs8094/Cleft/GSE59592/idat")
#list.files("panfs/panasas01/sscm/gs8094/Cleft/GSE59592/idat")

GSE62219 <- getGEOSuppFiles("GSE62219") #120/2 available
untar("/panfs/panasas01/sscm/gs8094/Cleft/GSE62219/GSE62219_RAW.tar", exdir = "/panfs/panasas01/sscm/gs8094/Cleft/GSE62219/idat")
list.files("panfs/panasas01/sscm/gs8094/Cleft/GSE62219/idat", pattern="idat")

GSE67444 <- getGEOSuppFiles("GSE67444") #140/2 available
untar("/panfs/panasas01/sscm/gs8094/Cleft/GSE67444/GSE67444_RAW.tar", exdir = "/panfs/panasas01/sscm/gs8094/Cleft/GSE67444/idat")
list.files("panfs/panasas01/sscm/gs8094/Cleft/GSE67444/idat")

#GSE36054 <- getGEOSuppFiles("GSE36054") #None available
#untar("panfs/panasas01/sscm/gs8094/Cleft/GSE36054/GSE36054_RAW.tar", exdir = "panfs/panasas01/sscm/gs8094/Cleft/GSE36054/idat")
#list.files("panfs/panasas01/sscm/gs8094/Cleft/GSE36054/idat")

#GSE36064 <- getGEOSuppFiles("GSE36064") #None available
#untar("panfs/panasas01/sscm/gs8094/Cleft/GSE36064/GSE36064_RAW.tar", exdir = "panfs/panasas01/sscm/gs8094/Cleft/GSE36064/idat")
#list.files("panfs/panasas01/sscm/gs8094/Cleft/GSE36064/idat")

GSE72556 <- getGEOSuppFiles("GSE72556") #196/2 available
untar("/panfs/panasas01/sscm/gs8094/Cleft/GSE72556/GSE72556_RAW.tar", exdir = "/panfs/panasas01/sscm/gs8094/Cleft/GSE72556/idat")
list.files("panfs/panasas01/sscm/gs8094/Cleft/GSE72556/idat")

#Untar all those that we managed to download
idatFilesGSE67444 <- list.files("/panfs/panasas01/sscm/gs8094/Cleft/GSE67444/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFilesGSE67444, gunzip, overwrite = TRUE)
idatFilesGSE62219 <- list.files("/panfs/panasas01/sscm/gs8094/Cleft/GSE62219/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFilesGSE62219, gunzip, overwrite = TRUE)
idatFilesGSE72556 <- list.files("/panfs/panasas01/sscm/gs8094/Cleft/GSE72556/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFilesGSE72556, gunzip, overwrite = TRUE)

#Get the pData (phenotype data) for the samples we've managed to download (226)
require(GEOquery)
geoMat <- getGEO("GSE67444",destdir="/panfs/panasas01/sscm/gs8094/Cleft/GSE67444") #70
pD.all <- pData(geoMat[[1]])
pD.GSE67444 <- pD.all[, c("title", "geo_accession", "characteristics_ch1", "characteristics_ch1.1","characteristics_ch1.2","characteristics_ch1.3", "characteristics_ch1.4","characteristics_ch1.5","characteristics_ch1.6")]

geoMat <- getGEO("GSE62219",destdir="/panfs/panasas01/sscm/gs8094/Cleft/GSE62219") #60
pD.all <- pData(geoMat[[1]])
pD.GSE62219 <- pD.all[, c("title", "geo_accession", "characteristics_ch1", "characteristics_ch1.1")]

geoMat <- getGEO("GSE72556",destdir="/panfs/panasas01/sscm/gs8094/Cleft/GSE72556") #96
pD.all <- pData(geoMat[[1]])
pD.GSE72556 <- pD.all[, c("title", "geo_accession", "characteristics_ch1", "characteristics_ch1.1","characteristics_ch1.2","characteristics_ch1.3", "characteristics_ch1.4","characteristics_ch1.5","characteristics_ch1.6")]

colnames(pD.GSE67444) <- c("title","geo_accession","ch_lead","mum_lead","sex","age_mnths","m_age_mnths","gestage_months","matsm")
colnames(pD.GSE62219) <- c("title","geo_accession","age_mnths","sex")
colnames(pD.GSE72556) <- c("title","geo_accession","adult_age_months","adult_bmi","adult_waist","sex","c_age_yrs","ch_bmi","ch_waist")

save(list=c("pD.GSE67444","pD.GSE62219","pD.GSE72556"),file="/panfs/panasas01/sscm/gs8094/Cleft/controls.phen.Robj")
load("/panfs/panasas01/sscm/gs8094/Cleft/controls.phen.Robj")

#Create the samplesheet
samplesheet.controls <- do.call("rbind",lapply(c(
  "/panfs/panasas01/sscm/gs8094/Cleft/GSE67444/idat",
  "/panfs/panasas01/sscm/gs8094/Cleft/GSE62219/idat"
  ),meffil.create.samplesheet))

samplesheet.controls <- merge(samplesheet.controls,pD.GSE67444,by.x="Sample_Name",by.y="geo_accession",all.x=T)
samplesheet.controls <- merge(samplesheet.controls,pD.GSE62219,by.x="Sample_Name",by.y="geo_accession",all.x=T)

#Clean up the samplesheet.controls data
library(stringr)
samplesheet.controls$Sex[samplesheet.controls$sex=="child gender: F"]<-"F"
samplesheet.controls$Sex[samplesheet.controls$sex=="child gender: M"]<-"M"
samplesheet.controls$Sex[samplesheet.controls$sex.x=="gender: F"]<-"F"
samplesheet.controls$Sex[samplesheet.controls$sex.x=="gender: M"]<-"M"

#c_age in weeks
samplesheet.controls$age_mnths.x <- as.numeric(unlist(lapply(sapply(as.character(samplesheet.controls$age_mnths.x),strsplit,split=": "),function(x) x[2])))*4
samplesheet.controls$age_mnths.y <- as.numeric(unlist(lapply(sapply(as.character(samplesheet.controls$age_mnths.y),strsplit,split=": "),function(x) x[2])))*4
#samplesheet.controls$c_age_yrs<- as.numeric(unlist(lapply(sapply(as.character(samplesheet.controls$c_age_yrs),strsplit,split=": "),function(x) x[2])))*52
samplesheet.controls$c_age <- rowSums(data.frame(samplesheet.controls$age_mnths.x ,samplesheet.controls$age_mnths.y 
#,samplesheet.controls$c_age_yrs
),na.rm=T)
samplesheet.controls$c_age <- ifelse(samplesheet.controls$c_age==0,NA,samplesheet.controls$c_age)

samplesheet.controls <- samplesheet.controls[,c("Sample_Name","Sex","Slide","sentrix_row","sentrix_col","Basename","c_age")]
samplesheet.controls$Study<-c(rep("GSE62219",60),rep("GSE67444",70))

samplesheet.controls$cleft_type<-"Control"
samplesheet.controls$CLO<-0
samplesheet.controls$CPO<-0
samplesheet.controls$CLP<-0
samplesheet.controls$CPO.CLO<-NA
samplesheet.controls$CPO.CLP<-NA
samplesheet.controls$CLO.CLP<-NA
samplesheet.controls$sample_type<-"Blood"
samplesheet.controls$participant <- NA

save(samplesheet.controls,file="samplesheet.controls.Robj")

############################################################
# Create full samplesheet (Cleft, Controls)   #
############################################################

#rbind together common column names
load("samplesheet.cleft.Robj")
load("samplesheet.controls.Robj")
Colnames<-intersect(colnames(samplesheet.controls),colnames(samplesheet.cleft))
samplesheet<-as.data.frame(rbind(samplesheet.controls[,Colnames],samplesheet.cleft[,Colnames]))
save(samplesheet,file="samplesheet.all.Robj")

#######################################################################################################################################

#############################################
#############################################
# Run Normalisation and QC with all samples #
#############################################
#############################################

setwd("Cleft/Controls2")
library(meffil)
load("samplesheet.all.Robj")


######################
# Perform sample QC  #
######################

qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=TRUE)
save(qc.objects,file="qc.objects.full.Robj")
load("qc.objects.full.Robj")
length(qc.objects) # check number of samples is correct
names(qc.objects) # sample names

#QC report

#Set some parameters for the QC of the raw data:
#beadnum.samples.threshold` = fraction of probes that failed the threshold of 3 beads.
#detectionp.samples.threshold` = fraction of probes that failed a detection.pvalue threshold of 0.01.
#beadnum.cpgs.threshold` = fraction of samples that failed the threshold of 3 beads.
#detectionp.cpgs.threshold` = fraction of samples that failed the detection.pvalue threshold of 0.01.
#sex.outlier.sd` = number of standard deviations to determine whether sample is sex outlier 
#snp.concordance.threshold` = concordance threshold to include snps to calculate sample concordance 
#sample.genotype.concordance.threshold` = concordance threshold to determine whether sample is outlier

qc.parameters <- meffil.qc.parameters(
    beadnum.samples.threshold             = 0.1,
    detectionp.samples.threshold          = 0.1,
    detectionp.cpgs.threshold             = 0.1, 
    beadnum.cpgs.threshold                = 0.1,
    sex.outlier.sd                        = 5,
    snp.concordance.threshold             = 0.95,
    sample.genotype.concordance.threshold = 0.8
)

# Summarise the QC analysis of the raw data (QC summary):

qc.summary <- meffil.qc.summary(qc.objects,parameters = qc.parameters)
save(qc.summary,file="qc.summary.Robj")
load("qc.summary.Robj")

# Make a QC report (html document in the working directory)
meffil.qc.report(qc.summary, output.file="qc.report.html")

# Remove bad samples
# Removing all outliers in:
#                            sample.name                            issue
#98                  200687170098_R03C01                     Sex mismatch
#260                 200687170135_R02C01                     Sex mismatch
#200687170143_R01C02 200687170143_R01C02       Methylated vs Unmethylated
#296                 200693590007_R02C02                     Sex mismatch
#309                 200693590022_R02C01                     Sex mismatch
#136                 200693590022_R04C02                     Sex mismatch
#141                 200693590023_R06C01                     Sex mismatch
#337                 200693590025_R03C01                     Sex mismatch
#369                 200693590033_R03C02                     Sex mismatch
#380                 200693590035_R06C01                     Sex mismatch
#387                 200693590037_R03C02                     Sex mismatch
#169                 200693590037_R04C02                     Sex mismatch
#391                 200693590038_R01C01                     Sex mismatch
#397                 200693590038_R06C01                     Sex mismatch
#399                 200693590056_R01C01                     Sex mismatch
#179                 200693590060_R01C02                     Sex mismatch
#413                 200693590060_R06C01                     Sex mismatch
#419                 200693590104_R05C01                     Sex mismatch
#429                 200693590105_R05C01                     Sex mismatch
#12472                        GSM1522959     Control probe (spec1.ratio1)
#12902                        GSM1522959      Control probe (spec1.ratio)
#GSM1522959                   GSM1522959                Detection p-value
#191                          GSM1522959                X-Y ratio outlier
#GSM1522960                   GSM1522960                Detection p-value
#12475                        GSM1522962     Control probe (spec1.ratio1)
#12905                        GSM1522962      Control probe (spec1.ratio)
#13335                        GSM1522962      Control probe (spec2.ratio)
#13765                        GSM1522962     Control probe (spec1.ratio2)
#GSM1522962                   GSM1522962                Detection p-value
#192                          GSM1522962                X-Y ratio outlier
#13766                        GSM1522963     Control probe (spec1.ratio2)
#GSM1522963                   GSM1522963                Detection p-value
#GSM1522968                   GSM1522968       Methylated vs Unmethylated
#GSM15229681                  GSM1522968                Detection p-value
#13799                        GSM1522996     Control probe (spec1.ratio2)
#GSM1522996                   GSM1522996                Detection p-value
#193                          GSM1522996                X-Y ratio outlier
#GSM1523010                   GSM1523010                Detection p-value
#GSM1523011                   GSM1523011                Detection p-value
#GSM1523013                   GSM1523013       Methylated vs Unmethylated
#12526                        GSM1523013     Control probe (spec1.ratio1)
#12956                        GSM1523013      Control probe (spec1.ratio)
#13386                        GSM1523013      Control probe (spec2.ratio)
#13816                        GSM1523013     Control probe (spec1.ratio2)
#GSM15230131                  GSM1523013                Detection p-value
#194                          GSM1523013                X-Y ratio outlier
#GSM1523016                   GSM1523016       Methylated vs Unmethylated
#58                           GSM1646850                     Sex mismatch
#10827                        GSM1646865 Control probe (spec2.G.34730329)
#76                           GSM1646891                     Sex mismatch



badsamples<-merge(samplesheet,qc.summary$bad.samples,by.x="Sample_Name",by.y="sample.name",all.x=F)
badsamples<-badsamples[which(badsamples$issue=="Sex mismatch"),c("Sample_Name","Sex","participant","sample_type","issue")]
badsamples<-badsamples[order(badsamples$participant),]

#200693590022_R04C02 and 200693590104_R05C01 should be removed as they are true sex mismatches

qc.summary$bad.samples$issue[which(qc.summary$bad.samples$sample.name %in% c("200693590022_R04C02","200693590104_R05C01"))]<-"Sex mismatch to be removed"

# except outliers based on control probes as there are often a few probes for each control type only. 
# For the control probe categories, meffil advice is to remove slides with dye bias Control probe (dye.bias) which are captured by the normalisation probes (normC + normG)/(normA + normT) or bisulfite conversion control probes (Control probe (bisulfite1) and Control probe (bisulfite2)).
# I left sex mismatches in as it's more likely the phenotype data on reported sex is wrong rather than the methylation data in cleft samples

outlier <- qc.summary$bad.samples
table(outlier$issue)
index <- outlier$issue %in% c("Control probe (dye.bias)", 
                              "Methylated vs Unmethylated",
                              "X-Y ratio outlier",
                              "Low bead numbers",
                              "Detection p-value",
                              "Sex mismatch to be removed",
                              "Genotype mismatch",
                              "Control probe (bisulfite1)",
                              "Control probe (bisulfite2)")

outlier <- outlier[index,]

#                            sample.name                      issue
#200687170143_R01C02 200687170143_R01C02 Methylated vs Unmethylated
#GSM1522959                   GSM1522959          Detection p-value
#191                          GSM1522959          X-Y ratio outlier
#GSM1522960                   GSM1522960          Detection p-value
#GSM1522962                   GSM1522962          Detection p-value
#192                          GSM1522962          X-Y ratio outlier
#GSM1522963                   GSM1522963          Detection p-value
#GSM1522968                   GSM1522968 Methylated vs Unmethylated
#GSM15229681                  GSM1522968          Detection p-value
#GSM1522996                   GSM1522996          Detection p-value
#193                          GSM1522996          X-Y ratio outlier
#GSM1523010                   GSM1523010          Detection p-value
#GSM1523011                   GSM1523011          Detection p-value
#GSM1523013                   GSM1523013 Methylated vs Unmethylated
#GSM15230131                  GSM1523013          Detection p-value
#194                          GSM1523013          X-Y ratio outlier
#GSM1523016                   GSM1523016 Methylated vs Unmethylated
#136                 200693590022_R04C02 Sex mismatch to be removed
#419                 200693590104_R05C01 Sex mismatch to be removed

length(qc.objects) #430
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects) #417
save(qc.objects,file="qc.objects.clean.Robj")

#Re-run QC summary and report on clean dataset
qc.summary <- meffil.qc.summary(qc.objects,parameters=qc.parameters)
save(qc.summary, file="qc.summary.clean.Robj")
meffil.qc.report(qc.summary, output.file="qc.report.clean.html")

#Extract raw beta matrix
load("qc.objects.clean.Robj")
betas.raw <- meffil.load.raw.data(qc.objects)
save(betas.raw,file="betas.raw.Robj")

#Extract detection P-values
load("qc.objects.clean.Robj")
detP <- meffil.load.detection.pvalues(qc.objects)
save(detP,file="detP.Robj")

##########################
# Normalise QC'd samples #
##########################

# set number of cores to use for parallelization
options(mc.cores=16)

# if necessary, load the qc.objects and qc.summary
#load("qc.objects.clean.Robj")
#length(qc.objects)
#load("qc.summary.clean.Robj")

# Estimate the number of PCs to use to adjust the methylation levels for technical effects.
# We can use 10-fold cross validation to estimate the residual variance after fitting n number of PCs.
# The residuals should consistently decrease with increasing numbers of PCs.
# The idea is to choose the number of PCs at which the residuals decreases abruptly.
# Sometimes there is more than one "elbow", in which case, meffil suggest choosing the one with the highest number of PCs.

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="pc.fit.pdf",height=6,width=6)

#set the number of PCs to use (from the plot, it looks like it's around 15, but difficult to tell and want to be stringent, so going with 25)

pc <- 25

# Perform functional normalisation
# Bad CpGs due to poor detection scores in the qc.summary are also removed.

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste0("norm.obj.pc",pc,".Robj"))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta,file=paste0("norm.beta.pc",pc,".Robj"))

# Generate normalisation report
# In the normalisation report, PCA will be performed on the control matrix and on the most variable probes.
# In addition, associations will be calculated between the first 10 PCs and the batch variables that were specified in the samplesheet.
# The association tests can be used to identify possible outliers. 
# For example, if Slide is one of your batch variables, it gives a p-value for each Slide rather than an overall p-value. 
# Poor slides can be identified and removed post-normalization.
# In the norm.parameters variable, you can set your batch variables eg. Slide, plate, and tissue etc., the number of extracted PCs in the control matrix, 
# the number of extracted PCs from the normalized betas, number of variable probes and the p-value threshold used for the association testing.
# It is important to code your batch variables as factors.

# Check:
str(norm.objects[[1]]$samplesheet)

# You change it by running a loop

for (i in 1:length(norm.objects)){
norm.objects[[i]]$samplesheet$Slide<-as.factor(norm.objects[[i]]$samplesheet$Slide)
norm.objects[[i]]$samplesheet$Sex<-as.factor(norm.objects[[i]]$samplesheet$Sex)
norm.objects[[i]]$samplesheet$sentrix_row<-as.factor(norm.objects[[i]]$samplesheet$sentrix_row)
norm.objects[[i]]$samplesheet$sentrix_col<-as.factor(norm.objects[[i]]$samplesheet$sentrix_col)
norm.objects[[i]]$samplesheet$sample_type<-as.factor(norm.objects[[i]]$samplesheet$sample_type)
norm.objects[[i]]$samplesheet$cleft_type<-as.factor(norm.objects[[i]]$samplesheet$cleft_type)
norm.objects[[i]]$samplesheet$Study<-as.factor(norm.objects[[i]]$samplesheet$Study)
}

# indicate batch variables
batch_var<-c("Slide","sentrix_row", "sentrix_col","sample_type","cleft_type","Study") 

#set normalisation parameters
norm.parameters <- meffil.normalization.parameters(
    norm.objects,
    variables=batch_var,
    control.pcs=1:10,
    batch.pcs=1:10,
    batch.threshold=0.01
)

# Run pcs, normalization summary and make normalisation report.

pcs <- meffil.methylation.pcs(norm.beta,probe.range=20000)
save(pcs,file="pcs.norm.beta.Robj")
norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs,parameters=norm.parameters)
meffil.normalization.report(norm.summary, output.file="normalization.report.html")



#Correct Sex to predicted Sex
samplesheet$Sex[which(samplesheet$Sample_Name %in% c("200693590007_R02C02",
"200693590025_R03C01",
"200693590023_R06C01",
"200693590060_R01C02",
"200687170135_R02C01",
"200693590105_R05C01",
"200693590038_R01C01",
"200693590038_R06C01",
"200693590033_R03C02",
"200693590037_R03C02" ,
"200687170098_R03C01",
"200693590037_R04C02",
"200693590022_R02C01",
"200693590060_R06C01",
"200693590035_R06C01",
"200693590056_R01C01" ))] <- ifelse(
samplesheet$Sex[which(samplesheet$Sample_Name %in% c("200693590007_R02C02",
"200693590025_R03C01",
"200693590023_R06C01",
"200693590060_R01C02",
"200687170135_R02C01",
"200693590105_R05C01",
"200693590038_R01C01",
"200693590038_R06C01",
"200693590033_R03C02",
"200693590037_R03C02" ,
"200687170098_R03C01",
"200693590037_R04C02",
"200693590022_R02C01",
"200693590060_R06C01",
"200693590035_R06C01",
"200693590056_R01C01" ))]=="M","F","M")


# Extract cell counts

cc<-t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc<-data.frame(IID=row.names(cc),cc)
write.table(cc,"houseman.cellcounts.txt",sep="\t",row.names=F,col.names=T,quote=F)
cc<-read.table("houseman.cellcounts.txt",header=T)

#Merge Cell Counts with samplesheet

samplesheet<-merge(samplesheet,cc,by.x="Sample_Name",by.y="IID",all.x=T)
save(samplesheet,file="samplesheet.cells.Robj")

#save different samplesheets for different tissues
samplesheet.blood <- droplevels(samplesheet[which(samplesheet$sample_type=="Blood"),])
samplesheet.lip <- droplevels(samplesheet[which(samplesheet$sample_type=="Lip"),])
samplesheet.palate <- droplevels(samplesheet[which(samplesheet$sample_type=="Palate"),])
samplesheet.eithertissue <- droplevels(samplesheet[which(samplesheet$sample_type=="Tissue"),])
samplesheet.tissue <- droplevels(rbind(samplesheet.lip,samplesheet.palate,samplesheet.eithertissue))

save(samplesheet.blood,file="samplesheet.cells.blood.Robj") #213
save(samplesheet.tissue,file="samplesheet.cells.tissue.Robj") #159
save(samplesheet.lip,file="samplesheet.cells.lip.Robj") #93
save(samplesheet.palate,file="samplesheet.cells.palate.Robj") #57

save(samplesheet,file="samplesheet.cells.Robj")

#Remove Outliers if necessary
rowIQR <- rowIQRs(norm.beta, na.rm = T)
row2575 <- rowQuantiles(norm.beta, probs = c(0.25, 0.75), na.rm = T)
maskL <- norm.beta< row2575[,1] - 3 * rowIQR 
maskU <- norm.beta > row2575[,2] + 3 * rowIQR 
norm.beta[maskL] <- NA
norm.beta[maskU] <- NA
save(norm.beta,file=paste0("norm.beta.pc",pc,".IQRtrimmed.Robj"))
