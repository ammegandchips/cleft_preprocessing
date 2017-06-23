#Code to clean and organise CC phenotype data for 450k analysis

setwd("Cleft/")
library(readstata13)
raw.phen <- read.dta13("/projects/Cleft_Collective/Gemma-Sharp/CC002-GS-20161102.dta")

# Names of columns in raw.phen:
# [1] "Filename"                         "FinalBarcode"                    
# [3] "matching_id"                      "cleft_type"                      
# [5] "SampleType"                       "cc21_db_wks_bld_smple"           
# [7] "cc21_db_wks_tis_smple"            "cc21_db_child_gender"            
# [9] "cc11_db_matage_at_delivery"       "ccpn11_blq_mthr_ethncty"         
#[11] "ccpn11_blq_matocc"                "ccpn11_blq_mated"                
#[13] "ccpn11_blq_matage_concep"         "ccpn11_blq_mat_smk_current_per_d"
#[15] "ccpn11_blq_mat_smk_current"       "ccpn11_blq_mat_smk_concep_per_da"
#[17] "ccpn11_blq_mat_smk_concep"        "ccpn11_blq_mat_parity"           
#[19] "ccpn11_blq_mat_gm_ethncty"        "ccpn11_blq_mat_gf_ethncty"  

phen <- data.frame(
	final_barcode=raw.phen$FinalBarcode,                   
	sentrix_id=raw.phen$Filename, 
	participant=raw.phen$matching_id,
	sample_type=raw.phen$SampleType,  
	Sex.tmp=raw.phen$cc21_db_child_gender,                     
	cleft_type= raw.phen$cleft_type,
	sample_weeks_blood= raw.phen$cc21_db_wks_bld_smple,
	sample_weeks_tissue= raw.phen$cc21_db_wks_tis_smple,
	m_age_delivery= raw.phen$cc11_db_matage_at_delivery,
	m_age_conception=raw.phen$ccpn11_blq_matage_concep,
	m_education=raw.phen$ccpn11_blq_mated,
	m_occupation=raw.phen$ccpn11_blq_matocc,
	m_smoking_current=raw.phen$ccpn11_blq_mat_smk_current,
	m_smoking_current_per_day=raw.phen$ccpn11_blq_mat_smk_current_per_d,
	m_smoking_conception=raw.phen$ccpn11_blq_mat_smk_concep,
	m_smoking_conception_per_day=raw.phen$ccpn11_blq_mat_smk_concep_per_d,    
	parity=raw.phen$ccpn11_blq_mat_parity
)


phen$m_age <- ifelse(phen$m_age_delivery=="NA",phen$m_age_conception,phen$m_age_delivery)
phen$c_age <- ifelse(phen$sample_weeks_blood=="NA",phen$sample_weeks_tissue,phen$sample_weeks_blood)
phen$Sex <- ifelse(phen$Sex.tmp==2,"F",NA)
phen$Sex <- ifelse(phen$Sex.tmp==1,"M",phen$Sex)
phen$Sex.tmp <- NULL

phen$cleft_type[which(phen$cleft_type==1)]<-"CLO"
phen$cleft_type[which(phen$cleft_type=="2")]<-"CPO"
phen$cleft_type[which(phen$cleft_type=="3")]<-"CLP"

phen$CLO<-ifelse(phen$cleft_type=="CLO",1,0)
phen$CLP<-ifelse(phen$cleft_type=="CLP",1,0)
phen$CPO<-ifelse(phen$cleft_type=="CPO",1,0)

phen$CPO.CLO[phen$CPO==1]<-1
phen$CPO.CLO[phen$CLO==1]<-0
phen$CPO.CLP[phen$CPO==1]<-1
phen$CPO.CLP[phen$CLP==1]<-0
phen$CLO.CLP<-ifelse(phen$CLO==1,1,NA)
phen$CLO.CLP<-ifelse(phen$CLP==1,0,phen$CLO.CLP)

phen<-droplevels(phen)

save(phen, file="phen.cleaned.Robj")


