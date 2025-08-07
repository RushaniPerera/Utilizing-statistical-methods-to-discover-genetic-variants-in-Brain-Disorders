library('data.table')
library(dplyr)

dise_path= "/work/long_lab/Rushani/IEDLMM/data/BPD/BDO_BARD_GRU_merge.match.b38.num.csv"
weight_file_path= "/work/long_lab/Rushani/IEDLMM/weight_files/8_exp_idp_weights.txt"


process_lifted_file<- function(dise_path){
  original_df=fread(file=dise_path, header=TRUE,fill = TRUE,sep=",",check.names = FALSE,select =c("CHR","LOC"))
  names(original_df)[names(original_df) == "CHR"] <- "CHROM"
  return(original_df)
}

read_weight_file <- function(weight_file_path){
  weight_df=read.table(file=weight_file_path,header = TRUE,fill = TRUE)
  weight_df <- subset(weight_df, select = c("CHROM", "POS","beta"))
  names(weight_df)[names(weight_df) == "POS"] <- "LOC" ## changing POS column name to LOC
  names(weight_df)[names(weight_df) == "beta"] <- "beta_1" ## changing beta column name to beta_1
  return(weight_df)
}

ordered_df=process_lifted_file(dise_path)
weight_df=read_weight_file(weight_file_path)
x=nrow(ordered_df)

get_new_weight_df=function(ordered_df,weight_df,x){
  df_new_weight=data.frame(LOC=ordered_df$LOC,CHROM=ordered_df$CHR,n_miss=numeric(x),beta=numeric(x),alpha=numeric(x),gamma=numeric(x)) #new dataframe to store weights
  df_new_weight$order <- seq_len(nrow(df_new_weight)) #add a new column so that we can get it n order of genotype file
  merged_df <- left_join(df_new_weight,weight_df, by = c("LOC","CHROM")) #left join into genotype locations
  merged_df <- merged_df[order(merged_df$order), ] #order according to order number
  merged_df<- merged_df[, !names(merged_df) %in% "order"] #remove order column
  merged_df['beta']=merged_df['beta_1'] # getting values for beta
  merged_df$beta[is.na(merged_df$beta) | is.null(merged_df$beta)] <- 0 #putting zero for beta if weight not available
  merged_df$n_miss <- 0
  merged_df$alpha <- 0
  merged_df$gamma <- 0
  merged_df <- merged_df[, !(names(merged_df) == "beta_1")]
}

final_weight=get_new_weight_df(ordered_df,weight_df,x)

write.table(final_weight, file = "/work/long_lab/Rushani/IEDLMM/weight_files/BPD/BPD_8_weights_adjusted.txt", sep = " ",row.names = FALSE, col.names = TRUE,quote = FALSE)


