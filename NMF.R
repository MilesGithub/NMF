
library(NMF)

gene_expression<-readRDS("D:/projects/TCGA/data/BRCA_mrna_df.rds")
gene_expression<-assay(gene_expression)
gene_expression<-data.frame(gene_expression)
colnames(gene_expression) <- gsub("\\.", "-", colnames(gene_expression))

clinical<-read.delim("D:/projects/TCGA/immune_subtypes/data/NIHMS958212-supplement-2.csv", sep=",")
clinical<-clinical[clinical$TCGA.Study == "BRCA",]

subtypes<-c("BRCA.LumA","BRCA.Basal","BRCA.Normal","BRCA.LumB","BRCA.Her2")

for (i in 1:length(subtypes)){
  subtype<-subtypes[4]
  print(subtype)
  clinical_sub<-clinical[clinical$TCGA.Subtype == subtype,]
  gene_expression_sub<-gene_expression[,colnames(gene_expression) %in% clinical_sub$TCGA.Participant.Barcode]

  set.seed(123)
  var_genes <- apply(gene_expression_sub, 1, var)
  high_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:2000]
  expr_matrix_filtered <- gene_expression_sub[high_var_genes, ]
  
  nmf_result <- nmf(expr_matrix_filtered, 2:5, nrun=10, seed=123456)

  pdf(paste0(subtype,"_NMF_rank_survey.pdf"), width=10, height=7)
  plot(nmf_result)
  dev.off()
  
  pdf(paste0(subtype, "_consensusmap.pdf"), width=15, height=10)
  consensusmap(nmf_result, labCol=NA, labRow=NA)
  dev.off()
  
}
