library(biomaRt)
library(tidyverse)
db <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)

# IDなどの情報
mart1 <-
  getBM(attributes = c("ensembl_gene_id", # Gene stable ID
                       "entrezgene_id", # NCBI gene (formerly Entrezgene) ID
                       "refseq_mrna", # RefSeq mRNA ID
                       "external_gene_name", # Gene name 
                       "description", # Gene description
                       "chromosome_name", # Chromosome/scaffold name
                       "start_position", # Gene start (bp)
                       "end_position", # Gene end (bp)
                       "strand", # Strand
                       "gene_biotype" # Gene type
  ), 
  mart = hg)

# GO情報
mart2 <-
  getBM(attributes = c("ensembl_gene_id", # Gene stable ID
                       "go_id", # GO term accession
                       "name_1006", # GO term name
                       "namespace_1003" # GO domain
  ), 
  mart = hg)


# Ensemble IDのベクトル作成
gene_id <- unique(mart1$ensembl_gene_id)

# Ensembe IDベクトルを使ってlapply
mart1_2 <- lapply(gene_id, function(x){
  
  # Ensemble IDのx番目のBioMart情報を取得
  tmp <- mart1[ mart1$ensembl_gene_id == x,]
  
  # もし複数行であれば情報をまとめる
  if(nrow(tmp) > 1){
    # 各列にuniqueな情報を\\で区切って文字をつなげる。
    tmp <- 
      apply(tmp, MARGIN = 2, FUN = function(y){
        # 列を抜き出して、重複除去、空セル/NAセルを除去
        y <- unique(y)
        y <- y[y != "" & !is.na(y)]
        
        # 残った要素を//で繋げて1つの要素にする。
        if(length(unique(y)) > 1){
          y <- paste(y, collapse = "//")
        }
        
        # 何も残っていない場合は空セルを入れておく
        if(length(y) == 0){
          y <- ""
        }
        return(y)
      })    
  }
  return(tmp)
})

# リストの各要素を縦に繋げてdata.frame作成
mart1_2 <- do.call(what = rbind, args = mart1_2)



# Gene Ontology情報を持っていない遺伝子は除く
tmp <- mart2[ mart2$go_id != "" & 
                mart2$namespace_1003 != "", ] # 404874 4

# GO IDとGO termを1列にまとめる
tmp$go_id_term <- paste(tmp$go_id, tmp$name_1006, sep = "_")
tmp$go_id <- NULL  # GO ID列は削除
tmp$name_1006 <- NULL  # GO term列は削除

# GO カテゴリを列に展開する
tmp_wide <- tidyr::pivot_wider(data = tmp,
                                   id_cols = ensembl_gene_id, # 残したい列
                                   names_from = namespace_1003, # 展開したいカテゴリを含む列
                                   values_from = go_id_term # 展開したい値を含む列
)
dim(tmp_wide) # [1] 26323     4

## data.frameに成形する。
# 空のdata.frame作成
GO <- data.frame(matrix(nrow = nrow(tmp_wide), ncol = 3),
                 row.names = tmp_wide$ensembl_gene_id)
colnames(GO) <- colnames(tmp_wide)[2:4]

# 1つのカテゴリ内に複数のGOがある場合は、//でつなげる。
for (i in 1:nrow(GO)) {
  GO[ i ,"molecular_function"] <-
    paste(unlist(tmp_wide$molecular_function[i]),collapse = "//")
  GO[ i ,"biological_process"] <-
    paste(unlist(tmp_wide$biological_process[i]),collapse = "//")
  GO[ i ,"cellular_component"] <-
    paste(unlist(tmp_wide$cellular_component[i]),collapse = "//")
}


# 行名をgene_id列として保存
GO$gene_id <- rownames(GO)
# left_joinで結合
gene_anno <- dplyr::left_join(x = mart1_2, y = GO, by = c("ensembl_gene_id"="gene_id"))

# 書き出し
write.csv(gene_anno, file = "gene_annotation.csv", row.names = F)