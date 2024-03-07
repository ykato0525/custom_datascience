setwd("/work/NGS/4.Reports")
# 発現ファイル読み込み
TPM <- read.table("genes_TPM.txt", sep="\t", row.names = 1)

# EnsembleID_GeneNameを_で区切り、EnsembleIDを新規列gene_idに保存
TPM$gene_id <- stringr::str_split(string = rownames(TPM), 
                                  pattern = "_", 
                                  simplify = TRUE)[,1]

# アノテーション情報読み込み
gene_anno <- read.csv("../gene_annotation.csv")

# 発現データにアノテーション情報をjoinする
TPM_anno <- dplyr::left_join(x = TPM, y = gene_anno, by = c("gene_id"="ensembl_gene_id"))

write.csv(TPM_anno, "genes_TPM_anno.csv")