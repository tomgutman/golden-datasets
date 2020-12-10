#New file to check SNV metrics for DREAM datasets results
#LIBRARIES

# sudo apt-get install libxml2 libxml2-dev

load.lib<-c("stringr", "viridis", "data.table","ggplot2","gridExtra","UpSetR",
            "caret","catspec","devtools","tidyr","plyr","dplyr",
            "rpart","rpart.plot","rattle","e1071","openssl","sqldf")



install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) {install.packages(lib,dependencies=TRUE)}
sapply(load.lib,library,character=TRUE)



 metrics_calculator_snv <- function(insilico_vcf, truth_insilico_vcf){
   
  truth_tab <- truth_insilico_vcf
  truth_tab_formated <- read.table(truth_tab)
  colnames(truth_tab_formated) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")
  truth_tab_formated$mut <- paste0(truth_tab_formated$"#CHROM","_",truth_tab_formated$POS,"_",truth_tab_formated$REF,"_",truth_tab_formated$ALT)

  insilico_tab <- insilico_vcf
  insilico_tab_formated <- read.table(insilico_tab, fill= TRUE, row.names=NULL)
  insilico_tab_formated <- insilico_tab_formated[,1:7] 
  colnames(insilico_tab_formated) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER")
  insilico_tab_formated <- insilico_tab_formated[insilico_tab_formated$FILTER == "PASS",]
  insilico_tab_formated$mut <- paste0(insilico_tab_formated$"#CHROM","_",insilico_tab_formated$POS,"_",insilico_tab_formated$REF,"_",insilico_tab_formated$ALT)
 
  options(scipen = 999)
  TP <- insilico_tab_formated[insilico_tab_formated$mut %in% truth_tab_formated$mut,]
  FP <- insilico_tab_formated[!insilico_tab_formated$mut %in% truth_tab_formated$mut,]
  FN <- truth_tab_formated[!truth_tab_formated$mut %in% insilico_tab_formated$mut,]
  Recall <- nrow(TP)/nrow(truth_tab_formated) 
  Precision <-nrow(TP)/(nrow(TP)+nrow(FP)) 
  F1_Score <- 2 * Precision * Recall / (Precision +Recall)
  
  metrics_X <- data.frame(round(Recall,3), round(Precision,3), round(F1_Score,3), nrow(TP), nrow(FP), nrow(FN))
  colnames(metrics_X) <- c("Recall","Precision","F1_Score","FP","TP","FN")
  metrics_X
  
  return(metrics_X)}


 #Set your own directory
 setwd('~/Downloads')
 truth_insilico_1_snv_vcf <-"./truth.snvs.synthetic.challenge.set1.vcf"
 truth_insilico_2_snv_vcf <-"./truth.snvs.synthetic.challenge.set2.vcf"
 truth_insilico_3_snv_vcf <-"./truth.snvs.synthetic.challenge.set3.vcf"
 
 
 insilico_1_snv_bsc <- "./insilico_1_snvs_bsc.vcf"
 insilico_1_result_snv_bsc <- metrics_calculator_snv(insilico_1_snv_bsc,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_bsc
 
 insilico_2_snv_bsc <- "./insilico_2_snvs_bsc.vcf"
 insilico_2_result_snv_bsc <- metrics_calculator_snv(insilico_2_snv_bsc,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_bsc

 insilico_3_snv_bsc <- "./insilico_3_snvs_bsc.vcf"
 insilico_3_result_snv_bsc <- metrics_calculator_snv(insilico_3_snv_bsc,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_bsc




 insilico_1_snv_hmf <- "./insilico_1_snvs_hmf.vcf"
 insilico_1_result_snv_hmf <- metrics_calculator_snv(insilico_1_snv_hmf,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_hmf
 
 insilico_2_snv_hmf <- "./insilico_2_snvs_hmf.vcf"
 insilico_2_result_snv_hmf <- metrics_calculator_snv(insilico_2_snv_hmf,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_hmf
 
 insilico_3_snv_hmf <- "./insilico_3_snvs_hmf.vcf"
 insilico_3_result_snv_hmf <- metrics_calculator_snv(insilico_3_snv_hmf,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_hmf
 
 
 
 insilico_1_snv_oicr <- "./insilico_1_snvs_oicr.vcf"
 insilico_1_result_snv_oicr <- metrics_calculator_snv(insilico_1_snv_oicr,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_oicr
 
 insilico_2_snv_oicr <- "./insilico_2_snvs_oicr.vcf"
 insilico_2_result_snv_oicr <- metrics_calculator_snv(insilico_2_snv_oicr,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_oicr
 
 insilico_3_snv_oicr <- "./insilico_3_snvs_oicr.vcf"
 insilico_3_result_snv_oicr <- metrics_calculator_snv(insilico_3_snv_oicr,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_oicr
 
 insilico_1_snv_charite <- "./insilico_1_snvs_charite.vcf"
 insilico_1_result_snv_charite <- metrics_calculator_snv(insilico_1_snv_charite,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_charite
 
 insilico_2_snv_charite <- "./insilico_2_snvs_charite.vcf"
 insilico_2_result_snv_charite <- metrics_calculator_snv(insilico_2_snv_charite,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_charite
 
 insilico_3_snv_charite <- "./insilico_3_snvs_charite.vcf"
 insilico_3_result_snv_charite <- metrics_calculator_snv(insilico_3_snv_charite,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_charite
 
 insilico_1_result_snv_bsc
 insilico_1_result_snv_charite
 insilico_1_result_snv_oicr
 insilico_1_result_snv_hmf
 
 insilico_2_result_snv_bsc
 insilico_2_result_snv_charite
 insilico_2_result_snv_oicr
 insilico_2_result_snv_hmf
 
 insilico_3_result_snv_bsc
 insilico_3_result_snv_charite
 insilico_3_result_snv_oicr
 insilico_3_result_snv_hmf
 
 
 
 
 #DEVEL
 
 insilico_1_snv_curie <- "./insilico_1_snvs_curie.vcf"
 insilico_1_result_snv_curie <- metrics_calculator_snv(insilico_1_snv_curie,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_curie
 
 insilico_2_snv_curie <- "./insilico_2_snvs_curie.vcf"
 insilico_2_result_snv_curie <- metrics_calculator_snv(insilico_2_snv_curie,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_curie
 
 insilico_3_snv_curie <- "./insilico_3_snvs_curie.vcf"
 insilico_3_result_snv_curie <- metrics_calculator_snv(insilico_3_snv_curie,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_curie

 
 insilico_1_result_snv_curie
 insilico_2_result_snv_curie
 insilico_3_result_snv_curie