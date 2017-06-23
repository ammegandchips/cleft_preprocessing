# cleft_preprocessing
Pre-processing 450k data from the Cleft Collective  

## Used in the paper
* 450k data were generated using blood and tissue (lip or palate) samples for 150 children from the Cleft Collective.
* Data were QC'd and normalised (functional normalization) by me using meffil and this script [Meffil_QC_450k_CC_No_Controls.R](Meffil_QC_450k_CC_No_Controls.R). These are the data that were used in the manuscript ['Distinct DNA methylation profiles in subtypes of orofacial cleft' Sharp et al. 2017 Clinical Epigenetics](http://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-017-0362-2)

### Preparing phenotype data, generating a smoking score, calculating epigenetic age
* Phenotype data were prepared and saved as phen.cleaned.Robj using [cleaning_pheno_450k_CC.R](cleaning_pheno_450k_CC.R)
* Maternal smoking status was estimated from cleft blood methylation data using [smoking_status_450k_CC.R](smoking_status_450k_CC.R). This is described in the paper, but basically I generated 4 scores
    1. using smoking-associated CpGs from PACE
    2. using smoking-associated CpGs from Reese
    3. using smoking-associated CpGs from PACE, adjusted for methylation PCs
    4. using smoking-associated CpGs from Reese, adjusted for methylation PCs
* I also generated the same 4 scores but using cleft tissue data instead of blood
* There wasn't really much between them, so I went for reese_blood_adj (i.e. number 4 in the list above)
* Epigenetic age and age acceleration was calculated using Horvath's website and this script [EpigeneticAge_450k_CC.R](EpigeneticAge_450k_CC.R)

## Trying to add a control group
* Using [Meffil_QC_450k_CC_Controls2.R](Meffil_QC_450k_CC_Controls2.R), I also attempted to normalise these samples alongside 'control' samples from GEO: GSE6744, GSE62219, GSE72556
* all the controls had some kind of difference to the CC samples:
    + GSE67444 Dried blood spots 6 to 60 months (confounded by batch, tissue, not sure on ethnicity)
    + GSE72556 age 3 to 5 Latino saliva (confounded by batch, tissue and ethnicity, age (too old))
    + GSE62219 PBLs age 3 to 60 months all girls (confounded by batch, tissue, not sure about ethnicity)
* I dropped GSE72556 before pre-processing because individuals in this study were all too old
* I described what happened when I tried to run an EWAS using these data (controls + cases) in the sup material of the paper. Basically, the batch effect was too large to get any interpretable results.
