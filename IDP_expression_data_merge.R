library(data.table)
library(dplyr)

expression_file_path="/work/long_lab/Rushani/IEDLMM/weight_files/ExpressionDirectedSNPsWeights_b38.txt"
img_file_path="/work/long_lab/Rushani/IEDLMM/weight_files/final_lifted_idp_data.txt"


read_idp_file <- function(img_file_path){
  idp_df=read.table(file=img_file_path, header=TRUE,fill = TRUE)
  idp_df$idp_beta=scale(idp_df$idp_beta)
  idp_df$idp_beta=abs(idp_df$idp_beta)
  return(idp_df)
}

read_expression_file <- function(expression_file_path){
  exprsn_df=read.table(expression_file_path, header=TRUE,fill = TRUE)
  setnames(exprsn_df, old = c('ps','chr','n_miss','alpha','beta','gamma'),
           new = c('POS','CHROM','n_miss','alpha','beta_e','gamma'))
  exprsn_df$beta_exp= (exprsn_df$alpha+exprsn_df$beta_e)
  exprsn_df$beta_exp=scale(exprsn_df$beta_exp)
 exprsn_df$beta_exp=abs(exprsn_df$beta_exp)
  return(exprsn_df)
}

idp_df=read_idp_file(img_file_path)

exprsn_df=read_expression_file(expression_file_path)

common_snps<-function(idp_df,exprsn_df){
  merged_df <-merge(exprsn_df,idp_df,by=intersect(names(idp_df),names(exprsn_df)),all.x = TRUE,all.y = TRUE)
  merged_df <- merged_df %>% rowwise() %>%
    mutate(beta = ifelse(is.na(idp_beta),beta_exp, ifelse(is.na(beta_exp), (8*idp_beta), sum(beta_exp+(8*idp_beta) ))))
  return(merged_df)
}

df=common_snps(idp_df,exprsn_df)

#Testing
subset(df,df$CHROM=='1' & df$POS=='100612278')
tail(df)
duplicates2 <- df[duplicated(df[c("CHROM", "POS")]) | duplicated(df[c("CHROM", "POS")], fromLast = TRUE), ]


write.table(df, file = "/work/long_lab/Rushani/IEDLMM/weight_files/8_exp_idp_weights.txt", sep = " ",row.names = FALSE, col.names = TRUE)


