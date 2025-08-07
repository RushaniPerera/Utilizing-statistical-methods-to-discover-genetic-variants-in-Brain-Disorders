library(dplyr)
library(readxl)

gene_ids_path= "BPD/cS2G_8_img_geneIDs.txt" ## Gene mapped file
disgen_path="BPD_C0005586_disease_gda_summary.xlsx" # This is downloadable from Disgenet website


get_gene_id <- function(gene_ids_path){
  gene_ids_df =read.table(file=gene_ids_path,header = TRUE,fill = TRUE,sep=",")
  gene_ids_df <- gene_ids_df[gene_ids_df$pvalue < 10^(-5), ]
  names(gene_ids_df)[names(gene_ids_df) == "gene_name"] <- "Gene" 
  return(gene_ids_df )
}

gene_id_df=get_gene_id(gene_ids_path)

get_disgen_data <- function(disgen_path){
  disgen_df =read_excel(disgen_path)
  return(disgen_df )
}
disgen_df= get_disgen_data(disgen_path)

#Genes identifies in BPD using EDLMM
no_of_genes_identified=length(unique(gene_id_df$Gene))
cat("no. of genes identified =",no_of_genes_identified)

#validating with Disgenet
bbb <-merge(gene_id_df,disgen_df[c("Score_gda","Gene")],by=intersect(names(gene_id_df),names(disgen_df))) ## intersect of expression and idp
no_of_validated_genes=length(unique(bbb$Gene))
cat("no. of genes validated =",no_of_validated_genes)


### Uniquely validated genes

IE_path= "BPD/cS2G_8_img_geneIDs.txt"
E_path="EDLMM_raw/BPD/cS2G_BPD_geneIDs.txt"


IE_df= get_gene_id(IE_path)
E_df=get_gene_id(E_path)

IE_unique_genes =unique(IE_df$Gene)
E_unqiue_genes=unique(E_df$Gene)

inter <- intersect(IE_unique_genes,E_unqiue_genes)
cat("no.of common genes identified from IEDLMM and EDLMM=",length(inter))

### IEDLMM ######
unique_identified_IEDLMM = setdiff(IE_unique_genes,inter)
cat("no. of uniquely identified genes from IEDLMM=",length(unique_identified_IEDLMM))

unique_validated_intersect_IEDLMM= intersect(unique_identified_IEDLMM,unique(disgen_df$Gene))
cat("no.uniquely validated genes from IEDLMM=",length(unique_validated_intersect_IEDLMM))

disgen_df_2= disgen_df[c("Gene","Score_gda")]
inter_bpd=disgen_df_2$Gene[disgen_df_2$Gene %in% unique_validated_intersect_IEDLMM]
disgen_df_3=disgen_df_2[disgen_df_2$Gene %in% inter_bpd,] ### Unique genes disgenet score


#### EDLMM #####
unique_identified_EDLMM= setdiff(E_unqiue_genes,inter)
cat("no. of uniquely identified genes from EDLMM=",length(unique_identified_EDLMM))

unique_validated_intersect_EDLMM= intersect(unique_identified_EDLMM,unique(disgen_df$Gene))
cat("no.uniquely validated genes from EDLMM=",length(unique_validated_intersect_EDLMM))


inter_bpd_edlmm=disgen_df_2$Gene[disgen_df_2$Gene %in% unique_validated_intersect_EDLMM]
disgen_df_4=disgen_df_2[disgen_df_2$Gene %in% inter_bpd_edlmm,]### Unique genes disgenet score - EDLMM


## Validation of the intersect
unique_validated_intersect= intersect(inter,unique(disgen_df$Gene))
cat("no.commonly validated genes=",length(unique_validated_intersect))










