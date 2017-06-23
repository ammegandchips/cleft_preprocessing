# Code to QC and normalise the Cleft Collective 450k data using meffil#
# October 12th 2016 #
# Gemma Sharp #

# Would not run interactively using qsub -I -l walltime=12:00:00 as could not access RDSF data via qsub #
# Have to run on login node until find better alternative #

#######################################################################################################################################

######################
# Create samplesheet #
######################
dir.create(file.path("Cleft/", "NoControls"), showWarnings = FALSE)
setwd("Cleft/NoControls/")
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
samplesheet.cleft$participant<-as.character(samplesheet.cleft$participant)
save(samplesheet.cleft,file="samplesheet.cleft.Robj")
samplesheet<-samplesheet.cleft

#######################################################################################################################################

#############################################
#############################################
# Run Normalisation and QC with all samples #
#############################################
#############################################

setwd("Cleft/NoControls/")
library(meffil)
load("samplesheet.cleft.Robj")
samplesheet<-samplesheet.cleft

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
#> qc.summary$bad.samples
#                            sample.name                              issue
#109                 200687170051_R05C01                  X-Y ratio outlier
#10                  200687170098_R03C01                       Sex mismatch
#130                 200687170135_R02C01                       Sex mismatch
#200687170143_R01C02 200687170143_R01C02         Methylated vs Unmethylated
#50                  200687170143_R01C02         Control probe (bisulfite1)
#350                 200687170143_R01C02         Control probe (bisulfite2)
#3350                200687170143_R01C02 Control probe (nonpoly.G.23663352)
#4250                200687170143_R01C02 Control probe (nonpoly.R.18773482)
#6050                200687170143_R01C02   Control probe (spec1.R.46779338)
#6650                200687170143_R01C02   Control probe (spec1.R.53740460)
#7850                200687170143_R01C02   Control probe (spec2.R.29662396)
#8150                200687170143_R01C02   Control probe (spec2.R.17661470)
#15                  200687170143_R01C02                  X-Y ratio outlier
#3671                200687170158_R06C01 Control probe (nonpoly.G.70645401)
#8796                200687170160_R06C02       Control probe (spec1.ratio1)
#9096                200687170160_R06C02        Control probe (spec1.ratio)
#9699                200693590007_R02C01       Control probe (spec1.ratio2)
#166                 200693590007_R02C02                       Sex mismatch
#179                 200693590022_R02C01                       Sex mismatch
#48                  200693590022_R04C02                       Sex mismatch
#53                  200693590023_R06C01                       Sex mismatch
#207                 200693590025_R03C01                       Sex mismatch
#239                 200693590033_R03C02                       Sex mismatch
#249                 200693590035_R05C01                  X-Y ratio outlier
#250                 200693590035_R06C01                       Sex mismatch
#257                 200693590037_R03C02                       Sex mismatch
#81                  200693590037_R04C02                       Sex mismatch
#261                 200693590038_R01C01                       Sex mismatch
#267                 200693590038_R06C01                       Sex mismatch
#269                 200693590056_R01C01                       Sex mismatch
#91                  200693590060_R01C02                       Sex mismatch
#283                 200693590060_R06C01                       Sex mismatch
#289                 200693590104_R05C01                       Sex mismatch
#299                 200693590105_R05C01                       Sex mismatch

badsamples<-merge(samplesheet,qc.summary$bad.samples,by.x="Sample_Name",by.y="sample.name",all.x=F)
badsamples<-badsamples[which(badsamples$issue=="Sex mismatch"),c("Sample_Name","Sex","participant","sample_type","issue")]
badsamples<-badsamples[order(badsamples$participant),]

#           Sample_Name Sex participant sample_type        issue
#18 200693590007_R02C02   F         018       Blood Sex mismatch
#22 200693590025_R03C01   F         018         Lip Sex mismatch
#21 200693590023_R06C01   M         041       Blood Sex mismatch 
#31 200693590060_R01C02   M         041      Palate Sex mismatch
#3  200687170135_R02C01   F         049       Blood Sex mismatch
#34 200693590105_R05C01   F         049         Lip Sex mismatch
#####20 200693590022_R04C02   M         052       Blood Sex mismatch
#28 200693590038_R01C01   F         082       Blood Sex mismatch
#29 200693590038_R06C01   F         082         Lip Sex mismatch
#####33 200693590104_R05C01   F         083       Blood Sex mismatch
#23 200693590033_R03C02   F         085       Blood Sex mismatch
#26 200693590037_R03C02   F         085         Lip Sex mismatch
#2  200687170098_R03C01   M         102      Palate Sex mismatch
#27 200693590037_R04C02   M         102       Blood Sex mismatch
#19 200693590022_R02C01   F         107      Palate Sex mismatch
#32 200693590060_R06C01   F         107       Blood Sex mismatch
#25 200693590035_R06C01   F         115      Palate Sex mismatch
#30 200693590056_R01C01   F         115       Blood Sex mismatch

# except outliers based on control probes as there are often a few probes for each control type only. 
# For the control probe categories, meffil advice is to remove slides with dye bias Control probe (dye.bias) which are captured by the normalisation probes (normC + normG)/(normA + normT) or bisulfite conversion control probes (Control probe (bisulfite1) and Control probe (bisulfite2)).
# I left sex mismatches in as it's more likely the phenotype data on reported sex is wrong rather than the methylation data in cleft samples, except for two samples where the tissue and blood were predicted different genders

qc.summary$bad.samples$issue[which(qc.summary$bad.samples$sample.name %in% c("200693590022_R04C02","200693590104_R05C01"))]<-"Sex mismatch to be removed"

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
#109                 200687170051_R05C01          X-Y ratio outlier #palate
#200687170143_R01C02 200687170143_R01C02 Methylated vs Unmethylated
#50                  200687170143_R01C02 Control probe (bisulfite1)
#350                 200687170143_R01C02 Control probe (bisulfite2) #lip
#15                  200687170143_R01C02          X-Y ratio outlier
#48                  200693590022_R04C02 Sex mismatch to be removed #blood
#249                 200693590035_R05C01          X-Y ratio outlier #blood
#289                 200693590104_R05C01 Sex mismatch to be removed #blood


length(qc.objects) #300
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects) #295
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

# Estimate the number of PCs to use to adjust the methylation levels for technical effects.
# We can use 10-fold cross validation to estimate the residual variance after fitting n number of PCs.
# The residuals should consistently decrease with increasing numbers of PCs.
# The idea is to choose the number of PCs at which the residuals decreases abruptly.
# Sometimes there is more than one "elbow", in which case, meffil suggest choosing the one with the highest number of PCs.

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="pc.fit.pdf",height=6,width=6)

#set the number of PCs to use (from the plot, it looks like it's around 15)

pc <- 15

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
}

# indicate batch variables
batch_var<-c("Slide","sentrix_row", "sentrix_col","sample_type","cleft_type") 

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

# Extract cell counts

cc<-t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc<-data.frame(IID=row.names(cc),cc)
write.table(cc,"houseman.cellcounts.txt",sep="\t",row.names=F,col.names=T,quote=F)

######################################
######################################
##  Add extra data to samplesheet   ##
######################################
######################################

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

#Merge Cell Counts with samplesheet
cc<-read.table("houseman.cellcounts.txt",header=T)
samplesheet<-merge(samplesheet,cc,by.x="Sample_Name",by.y="IID",all.x=T)
save(samplesheet,file="samplesheet.cleft.Robj")

#Add smoking scores
load("Smoking_Status_Score.Rdata")
Scores<-All
samplesheet<-samplesheet[order(samplesheet$participant,samplesheet$sample_type),]
samplesheet<-merge(samplesheet,Scores,by="participant",all=T)
samplesheet$m_smoking_conception[which(samplesheet$m_smoking_conception==2)]<-0
samplesheet$m_smoking_current[which(samplesheet$m_smoking_current==2)]<-0

samplesheet$smoking_missing<-0
samplesheet$smoking_missing[is.na(samplesheet$m_smoking_conception)]<-"missing"
samplesheet$smoking_missing[which(samplesheet$m_smoking_conception==1)]<-"smoker"
samplesheet$smoking_missing[which(samplesheet$m_smoking_conception==0)]<-"non-smoker"

#save different samplesheets for different tissues
samplesheet.blood <- droplevels(samplesheet[which(samplesheet$sample_type=="Blood"),])
samplesheet.lip <- droplevels(samplesheet[which(samplesheet$sample_type=="Lip"),])
samplesheet.palate <- droplevels(samplesheet[which(samplesheet$sample_type=="Palate"),])
samplesheet.eithertissue <- droplevels(samplesheet[which(samplesheet$sample_type=="Tissue"),])
samplesheet.tissue <- droplevels(rbind(samplesheet.lip,samplesheet.palate,samplesheet.eithertissue))

save(samplesheet.blood,file="samplesheet.blood.Robj") #213
save(samplesheet.tissue,file="samplesheet.tissue.Robj") #159
save(samplesheet.lip,file="samplesheet.lip.Robj") #93
save(samplesheet.palate,file="samplesheet.palate.Robj") #57

samplesheet->samplesheet.cleft
save(samplesheet.cleft,file="samplesheet.cleft.Robj")

#remove outliers using IQR3 (Tukey) approach
rowIQR <- rowIQRs(norm.beta, na.rm = T)
row2575 <- rowQuantiles(norm.beta, probs = c(0.25, 0.75), na.rm = T)
maskL <- norm.beta< row2575[,1] - 3 * rowIQR 
maskU <- norm.beta > row2575[,2] + 3 * rowIQR 
norm.beta[maskL] <- NA
norm.beta[maskU] <- NA
save(norm.beta,file=paste0("norm.beta.pc",pc,".IQRtrimmed.Robj"))
