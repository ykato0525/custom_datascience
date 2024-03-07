#!/bin/bash


# 開始時間を表示
echo starting at
date

# 解析環境に移動する
cd /work/NGS

# エラー回避のための上限を挙げておく
ulimit -n 16000

# 変数の設定
STAR_REF=/work/NGS/Reference/STAR_index
RSEM_REF=/work/NGS/Reference/RSEM_Reference

# サンプル名とその他の情報を分ける区切り文字
DELIM=_

# BAM indexファイルを作るかどうか。YES以外であれば作られない。
CREATE_BAI=YES

# ディレクトリ構造を整理する
mkdir {0.QC,0.QC/fastqc,1.fastq,1.fastq/raw,1.fastq/after_fastp,2.STAR,2.STAR/BAM,3.RSEM,4.Reports}



# ディレクトリのパスを変数に保存する
WD=$(pwd)
QC=${WD}/0.QC
FQ_RAW=${WD}/1.fastq/raw
FQ=${WD}/1.fastq/after_fastp
STAR=${WD}/2.STAR
RSEM=${WD}/3.RSEM
REPORTS=${WD}/4.Reports

mv *.fastq.gz ./1.fastq/raw

# ファイルの情報を取得する
FILENAME=$(ls ${FQ_RAW}/*{fq.gz,fastq.gz} | xargs basename -a | sed "s/\(.*\)${DELIM}.*/\1/" | uniq)
echo ----- Analyse files: -----
echo $FILENAME
READ1=$(ls ${FQ_RAW}/*{fq.gz,fastq.gz} | xargs basename -a | sed "s/\(.*\)${DELIM}\(.*\)/\2/" | sort | uniq | awk 'NR==1')
READ2=$(ls ${FQ_RAW}/*{fq.gz,fastq.gz}| xargs basename -a | sed "s/\(.*\)${DELIM}\(.*\)/\2/" | sort | uniq | awk 'NR==2')
echo ----- Extension: -----
echo $READ1 $READ2


# fastqcを行う
echo "start fastqc"
## fastqc for raw fastq ##
fastqc \
-t 16 \
--nogroup \
--quiet \
-o ${QC}/fastqc \
${FQ_RAW}/*{fq.gz,fastq.gz}

##### fastp #####
for i in $FILENAME
do
fastp -i ${FQ_RAW}/${i}_${READ1} \
-I ${FQ_RAW}/${i}_${READ2} \
-o ${FQ}/${i}_trim_${READ1} \
-O ${FQ}/${i}_trim_${READ2} \
-w 16 \
--detect_adapter_for_pe \
--length_required 25 \
-j ${FQ}/${i}_fastp.json \
-h ${FQ}/${i}_fastp.html
done

echo "fastp done" 

# STAR 
# bam indexファイルを作るかどうかでSTARの--outSAMtype引数を変更する。
if [ $CREATE_BAI = YES ]; then \
OUTSAMTYPE="BAM SortedByCoordinate"
mkdir ${REPORTS}/bam
else
OUTSAMTYPE="None"
fi

for i in $FILENAME
do
STAR --genomeDir ${STAR_REF} \
--runThreadN 16 \
--outFileNamePrefix ${STAR}/${i}_ \
--quantMode TranscriptomeSAM \
--outSAMtype ${OUTSAMTYPE} \
--readFilesCommand zcat \
--readFilesIn ${FQ}/${i}_trim_${READ1} ${FQ}/${i}_trim_${READ2}
done

if [ $CREATE_BAI = YES ]; then \
for i in $FILENAME ; do \
mv ${STAR}/${i}_Aligned.sortedByCoord.out.bam ${REPORTS}/bam/${i}.bam \
; samtools index -@ 16 ${REPORTS}/bam/${i}.bam ${REPORTS}/bam/${i}.bam.bai; done
fi

# RSEM
for i in $FILENAME
do
rsem-calculate-expression -p 16 \
--paired-end \
--alignments \
--append-names \
--estimate-rspd \
--no-bam-output \
--quiet \
${STAR}/${i}_Aligned.toTranscriptome.out.bam \
${RSEM_REF} \
${RSEM}/${i}
done

# multiqc for STAR
multiqc ${STAR}/ -o ${QC}/STAR


# merge all file
cd ${RSEM}

cut -f 1 $( echo "$FILENAME" | awk 'NR==1').genes.results > genes.join.txt

for i in $FILENAME
do
join -a 1 -a 2 -j 1 -e NA \
<(head -n 1 genes.join.txt && tail -n +2 genes.join.txt | sort -k1,1) \
<(head -n 1 ${i}.genes.results | sed -e "s/\t/\t${i}_/g" && tail -n +2 ${i}.genes.results | sort -k1,1)\
| tr " " "\t" > genes.join2.txt \
; cp genes.join2.txt genes.join.txt \
; rm genes.join2.txt
done

cut -f 1 $( echo "$FILENAME" | awk 'NR==1').isoforms.results > isoforms.join.txt

for i in $FILENAME
do
join -a 1 -a 2 -j 1 -e NA \
<(head -n 1 isoforms.join.txt && tail -n +2 isoforms.join.txt | sort -k1,1) \
<(head -n 1 ${i}.isoforms.results | sed -e "s/\t/\t${i}_/g" && tail -n +2 ${i}.isoforms.results | sort -k1,1) \
| tr " " "\t" > isoforms.join2.txt \
; cp isoforms.join2.txt isoforms.join.txt \
; rm isoforms.join2.txt
done

## TPM列, expected count列だけのファイルを作成
# gene expression TPM
cut -f 1,$(head -1 genes.join.txt | xargs -n1 | sed -n "/TPM/=" | paste -s | tr "\t" ",") genes.join.txt | sed "1s/_TPM//g" > genes_TPM.txt
# gene expression expected count
cut -f 1,$(head -1  genes.join.txt| xargs -n1 | sed -n "/expected_count/=" | paste -s | tr "\t" ",") genes.join.txt  | sed "1s/_expected_count//g" > genes_expected_count.txt
# transcript expression TPM
cut -f 1,$(head -1  isoforms.join.txt | xargs -n1 | sed -n "/TPM/=" | paste -s | tr "\t" ",") isoforms.join.txt | sed "1s/_TPM//g" > isoforms_TPM.txt
# transcript expression expected count
cut -f 1,$(head -1 isoforms.join.txt | xargs -n1 | sed -n "/expected_count/=" | paste -s | tr "\t" ",") isoforms.join.txt  | sed "1s/_expected_count//g" > isoforms_expected_count.txt

## 成型したファイルを4.Reportsフォルダに移動
mv genes_TPM.txt genes_expected_count.txt isoforms_TPM.txt isoforms_expected_count.txt ${REPORTS}/.

cd ${WD}

Rscript join_annotation.R

echo ending at
date