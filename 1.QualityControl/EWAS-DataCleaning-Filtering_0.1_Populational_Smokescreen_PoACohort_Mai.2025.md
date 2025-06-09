# EWAS - Clinical Controls (Populational 2) - Data cleaning and filtering

author: "Iago"

date: "2025-02-03"

This script is dedicated to filter individuals from Populational-Smokescreen and also add the Principal Component data to the SPSS database.

#### purl this file to Rscript

```{r}

scriptsFolder = "/media/iago/Backup/Ubuntu/Documentos/1.projetos_e_resultados/EWAS-proj/3.EWAS-PoA-ADHD_Pipeline/1.Populational/Mai.2025/4.Subset-Smokescreen/5.Populational_Smokescreen_PoACohort_Rscripts/1.QualityControl/"

Rmd = file.path(scriptsFolder,"EWAS-DataCleaning-Filtering_0.1_Populational-Smokescreen-PoACohort_Mai.2025.Rmd")

Rscript = file.path(scriptsFolder,"EWAS-DataCleaning-Filtering_0.1_Populational-Smokescreen-PoACohort_Mai.2025.R")
 knitr::purl(Rmd, output = Rscript)

```

``` {.bash .text}
PCA_Smokescreen_Data="/scratch/iago/EWAS/3.EWAS-PoA-ADHD/1.Populacional/Mai.2025/4.Subset-Smokescreen/4.Original_MetaData/GWAS-Smokescreen_PoACohort_Mar.2025_PCA-Covariates.eigenvec.header"
```

# 

```{r}

source("/scratch/iago/EWAS/3.EWAS-PoA-ADHD/source_functions/source_functions.R")

need_packages <-c("tibble", "minfi", "readxl", "stringr", "dplyr","haven")#' install by running the function
install_cran_pkgs(pkgs = need_packages)
install_bioconductor_pkgs(pkgs = need_packages)

#' Check if all you need is installed
check_installed(pkgs = need_packages)

#load all packages, if needed
load_pkgs <- c("minfi", "readxl", "stringr", "dplyr","haven")
lapply(load_pkgs, require, character.only = TRUE )
```

# 1. Joining correspondence sheet with clinical IDS

```{r}

corr = read.table("/scratch/iago/GWAS-PoACohort/LAGC-Pipeline/GWAS-Smokescreen/Mar.2025/1.RawData/Smokescreen_anotacao_fam.txt", sep=" ", header =F)


pheno = read.csv("/scratch/iago/EWAS/3.EWAS-PoA-ADHD/1.Populacional/Mai.2025/4.Subset-Smokescreen/6.Phenotype_files/PhenoFile-Populational_Smokescreen_PoACohort_Mai.2025.csv")


pheno$IID = as.numeric(str_extract(pheno$smoke_FID, "(?<=\\*)\\d+"))

pheno <- dplyr::left_join(pheno, corr, by = c("IID" = "V2"))

pheno2 = pheno %>% dplyr::select(c(
  "barcode", "smoke_FID", "Sentrix_ID", "Sentrix_Position",
  "EWAS", "TDAH", "Sex", "Age",
  "Smoking", "ASRS_DESmedia", "ASRS_HIPIMPmedia", "ASRS_HIPmedia",
  "ASRS_IMPmedia", "ASRS_TDAH", "IID","V7"
)
) %>% rename(IID_CNB = V7)

write.csv(pheno2, "/scratch/iago/EWAS/3.EWAS-PoA-ADHD/1.Populacional/Mai.2025/4.Subset-Smokescreen/4.Original_MetaData/PhenoFile-CORRESPONDENCE-Populational_Smokescreen_PoACohort_Mai.2025.csv", row.names = F, quote = F)


```

# 2. Selecting individuals and saving phenoFile

```{r}

## 0.1 Directory variables -----------------------------------------------------
baseFolder="/scratch/iago/EWAS/3.EWAS-PoA-ADHD/1.Populacional/Mai.2025/4.Subset-Smokescreen/"

inputFolder <- "/scratch/physiogenlab/EWASdata/input_data/EWAS_PoA_cohort_population/all_idat_files_Jan.2024/"

MetaData_inputFolder=file.path(baseFolder,"4.Original_MetaData/")

outputFolder= file.path(baseFolder,"6.Phenotype_files/")
## 0.1 Input variables  ---------------------------------------------
banco <- read_sav(paste0(MetaData_inputFolder,"Banco_Geral_TDAH-CONTROLES_EWAS_23.03.24.sav"))
date <- "Mai.2025"

# - basenameForFiles
basenameForFiles <- "Populational_Smokescreen_PoACohort"

#Samplesheet original : Tabela contendo a relação FID IID --> SentrixID SentrixPosition

samplesheet <- as.data.frame(read_excel(paste0(inputFolder,"Populacional_sample_RS_Brazil_covar.xlsx")))
samplesheet <-  samplesheet %>% select(SampleID,ID_populacional)


# PCA - Plink genetic PCs

pca = read.table(file.path(baseFolder, "4.Original_MetaData/GWAS-Smokescreen_PoACohort_Mar.2025_PCA-Covariates.eigenvec.header"), 
                 sep = " ", header = T)


correspondence = read.csv(file.path(baseFolder, "4.Original_MetaData/PhenoFile-CORRESPONDENCE-Populational_Smokescreen_PoACohort_Mai.2025.csv"))
                                      
                                      
## 0.3 Output variables --------------------------------------------------------
populational_Smokescreen_phenoFile.file <- paste0(outputFolder,"PhenoFile-",basenameForFiles,"_",date,".csv")

```

```{r}

#' Encontrar 5 individuos casos (aleatórios) e copiar seus idat files para 
#' a pasta:
#' "/home/physiogenlab/EWASdata/input_data/EWAS_PoA_cohort_clinical/3.Controles_Populational-1"


banco_Smokescreen <- banco %>% filter(EWAS ==2, batch_EWAS_array == 1) %>% 
  select(  c("SampleID",
             "smoke_FID",
             "SentrixID", 
             "SentrixPosition", 
             "EWAS", 
             "TDAH",
             "sexo_nov23", 
             "idade_baseline_nov23",
             "usonicot_nov23",
             "ASRSfromsnap_DESmedia_corr",
             "ASRSfromsnap_HIPIMPmedia_corr",
             "ASRSfromsnap_HIPmedia_corr",
             "ASRSfromsnap_IMPmedia_corr",
             "ASRSfromsnap_TDAH",
)) %>%
  rename(barcode = SampleID, 
         Sentrix_ID = SentrixID, 
         Sentrix_Position = SentrixPosition, 
         Sex = sexo_nov23, 
         Smoking = usonicot_nov23, 
         ASRS_DESmedia = ASRSfromsnap_DESmedia_corr,
         ASRS_HIPmedia = ASRSfromsnap_HIPmedia_corr,
         ASRS_HIPIMPmedia = ASRSfromsnap_HIPIMPmedia_corr,
         ASRS_IMPmedia = ASRSfromsnap_IMPmedia_corr,
         ASRS_TDAH = ASRSfromsnap_TDAH,
         Age = idade_baseline_nov23)

corre = correspondence %>% dplyr::select(barcode,IID_CNB)
banco_Smokescreen2 = dplyr::left_join(banco_Smokescreen, corre, by = "barcode")

pca_5PCs = pca %>% dplyr::select(IID,PC1,PC2,PC3,PC4,PC5)

banco_Smokescreen_PCs = dplyr::left_join(banco_Smokescreen2, pca_5PCs, by = c("IID_CNB" = "IID")) %>% 
  as.data.frame() %>% 
  dplyr::select("barcode", 
                "IID_CNB", 
                "Sentrix_ID", 
                "Sentrix_Position",
                "EWAS", 
                "TDAH", 
                "Sex", 
                "Age",
                "Smoking", 
                "ASRS_DESmedia", 
                "ASRS_HIPIMPmedia", 
                "ASRS_HIPmedia",
                "ASRS_IMPmedia", 
                "ASRS_TDAH", 
                "PC1", "PC2", "PC3", "PC4", "PC5" ) %>% 
  rename(IID=IID_CNB)
```

```{r}


write.csv(banco_Smokescreen_PCs, populational_Smokescreen_phenoFile.file, row.names = F)
```

Generate list of idat files to trasnfer to the new input folder:

```{r}

idat.listFiles <-  as.list(populational1$barcode)

writeLines(unlist(idat.listFiles), paste0(inputFolder,"idatFiles.list_",basenameForFiles,"_",date,".txt"))

```

------------------------------------------------------------------------

### 2. Move all the idat files to the new input folder

```{bash}
# 0.1 input variables ------------------------------------------------------------------------
idatFileList="/scratch/physiogenlab/EWASdata/input_data/EWAS_PoA_cohort_population/all_idat_files_Jan.2024/idatFiles.list_Populational_Smokescreen_PoACohort_Mai.2025.txt"

source_dir="/scratch/physiogenlab/EWASdata/input_data/EWAS_PoA_cohort_population/all_idat_files_Jan.2024/"

# 0.2 Output variables ------------------------------------------------------------------------
symLink_dir="/scratch/physiogenlab/EWASdata/input_data/EWAS_PoA_cohort_population/IdatFiles-Subset-Smokescreen_Mai.2025"

#-----------------------------------------------------------------------

while IFS= read -r line
do
  #echo "${line}_*.idat"
  # - Create a hard link of each .idat file present in the file list.
  ln "$source_dir/${line}"_*.idat $symLink_dir
done < "$idatFileList"

echo "Done."
```
