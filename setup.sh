#!/bin/bash

# NGSフォルダを作るようのbashファイルになります。
mkdir NGS/Refernce
# ReferenceフォルダにNGS関連のファイルを放り込む

# scriptたちを移動する
mv chip.sh NGS
mv join_annotation.R NGS
mv make_annotation.R NGS


cd NGS

# Referenceフォルダにリファレンス配列をインストール
cd Refernce
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.111.gtf.gz

# STARのセットアップ
STAR --runMode genomeGenerate \
--genomeDir /work/NGS/Reference/STAR_index \
--runThreadN 16 \
--genomeFastaFiles /work/NGS/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /work/NGS/Reference/Homo_sapiens.GRCh38.111.gtf

# RSEMのセットアップ
rsem-prepare-reference \
--num-threads 16 \
--gt/work/NGS ~/Reference/Homo_sapiens.GRCh381109.gtf /work/NGS
~/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa/work/NGS\
~/Reference/RSEM_Reference
