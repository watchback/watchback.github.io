---
# multilingual page pair id, this must pair with translations of this page. (This name must be unique)
lng_pair: id_RNA-seq
title: "[Bioinformatics] De Novo RNA-seq Analysis pipeline"

# post specific
# if not specified, .name will be used from _data/owner/[language].yml
author: HANJIN KIM
# multiple category is not supported
category: HANJIN KIM
# multiple tag entries are possible
tags: [Bioinformatics, RNA-seq, Pipeline, De Novo]
# thumbnail image for post
img: ":RNA-seq.jpg"
# disable comments on this page
#comments_disable: true

# publish date
date: 2023-12-28 10:24:30 +0900

# seo
# if not specified, date will be used.
#meta_modify_date: 2022-01-01 10:04:30 +0900
# check the meta_common_description in _data/owner/[language].yml
#meta_description: ""

# optional
# please use the "image_viewer_on" below to enable image viewer for individual pages or posts (_posts/ or [language]/_posts folders).
# image viewer can be enabled or disabled for all posts using the "image_viewer_posts: true" setting in _data/conf/main.yml.
#image_viewer_on: true
# please use the "image_lazy_loader_on" below to enable image lazy loader for individual pages or posts (_posts/ or [language]/_posts folders).
# image lazy loader can be enabled or disabled for all posts using the "image_lazy_loader_posts: true" setting in _data/conf/main.yml.
#image_lazy_loader_on: true
# exclude from on site search
#on_site_search_exclude: true
# exclude from search engines
#search_engine_exclude: true
# to disable this page, simply set published: false or delete this file
#published: false
---
<!-- outline-start -->

# **0) Raw data 및 환경설정**
#### 필요한 software  

- Trimmomatic  

- Trinity  

- Bowtie2  

- Samtools  

- TransDecoder  

- BLAST  

- CD-hit  

- RSEM  

- pb-assembly   

- pbbam   

- BUSCO   

- pilon   

- R  

- R packages (edgeR, ctc, Biobase, ape, gplots 등)

- Python  

bioconda 설치 or 직접 설치, PATH 설정 필요  
  
Bioconda를 이용해 설치하면 간편하나 가끔 Bioconda를 통해 설치하여도 동작하지 않는 경우 직접 설치 하는 방법도 알아야 함
  
  
---
#### RAW DATA
  
FASTQ 파일에는 sequence 서열 뿐만 아니라 quality 정보도 포함하고 있는데  

첫번째 줄은 “@” 로 시작하며 Sequence ID  

두번째 줄은 Sequence 서열  

세번째 줄은 “+” 하나만 있거나, 또는 그 뒤에 첫번째줄의 Sequence ID 의 반복  

네번째 줄은 각 서열의 quality 를 나타내는 기호 (ASCII code)로 구성  
  

파일이 압축되어 있는 경우 압축 해제  

File:**Forward Sequence.gz, Reverse Sequence.gz**  
Tool:**Linux command**  

#### 커맨드

tar.gz 

    tar -zxvf rawdata.fastq_1.tar.gz

gz 

    gzip -d rawdata.fastq_1.gz

#### 결과
**rawdata_1.fastq**  
**rawdata_2.fastq**  


# **1) Preprocessing(Quality Control)**
  
Illumina Sequencing 같은 경우는 Adapter 제거가 필요하나 PacBio Sequencing은 Adapter를 제거한 데이터를 주기 때문에 하지 않아도 됨  

File: **Forward Sequence.fastq, Reverse Sequence.fastq**  
Tool: **Trimmomatic**  

#### 커맨드 (.sh 파일로 만들어 실행)

    java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq input_reverse.fq output_forward_paired.fq output_forward_unpaired.fq output_reverse_paired.fq output_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#### 옵션
- -threads: thread 수를 지정   
  
- ILLUMINACLIP: fastaWithAdaptersEtc : seed mismatches : palindrome clip threshold : simple clip threshold  
  
fastaWithAdaptersEtc == adaptor 파일, Trimmomatic 프로그램 내에 adaptor 디렉토리에 존재  
  
seed mismatches == 최초 16 bases를 seed로 놓고 이를 full match가 허용하는만큼 확장한다. full match에 허용할 mismatch의 최소값(2)  
  
palindrome clip threshold == paried-ended data의 경우 score 값(30) /약 50 bases  
simple clip threshold == single-ended data의 경우 score 값(10) /약 17 bases  
  
- SLIDINGWINDOW:[num1]:[num2]: 주어진 window(similar to mer) 상수값[num1]만큼 sequence를 sliding하며 window falls [num1]내의 average quality가 주어진 값[num2]보다 낮을 경우 제거   
  
- LEADING (LEADING:3): a read의 앞쪽이 주어진 threshold quality보다 낮은 경우 제거  
  
- TRAILING (TRAILING:3): a read의 뒤쪽이 주어진 threshold quality보다 낮은 경우 제거  
  
- CROP : 명시된 길이만큼 read의 뒤쪽부분을 제거  
  
- HEADCROP : 명시된 길이만큼 read의 앞쪽부분을 제거  
  
- MINLEN:(num) read trimming 과정중에 (num)bp 미만의 read 는 버리라는 명령인데, sequencing data 가 101bp 인 경우에는 MINLEN:36, 151bp 인 경우에는 MINLEN:50을 줌  
  
- TOPHRED33 : quality scores를 Phread-33으로 변경  
    
- TOPHRED64 : quality scores를 Phread-64로 변경  
   

#### 결과
**output_forward_paired.fq**  
**output_reverse_paired.fq**  
**output_forward_unpaired.fq**  
**output_reverse_unpaired.fq**  

# **2) Concatenation of samples data**

De novo RNA-seq 분석을 하게 된다면 일반적으로 여러 samples(주로 multiple tissues)을 가지고 시작
genome 이 없다보니 genome 과 유사한 역할을 할 수 있는 transcriptome assembly 제작 필요  
이때 최대한 많은 transcriptome 정보가 담긴 transcriptome assembly 를 만들기 위해 여러 sample의 RNA-seq data를 한개의 파일로 만들어 줌
forward끼리, reverse끼리 concat  


File: **Forward Sequence.fastq, Reverse Sequence.fastq**  
Tool: **Linux command**  

#### 커맨드

    cat Forward Sequence1.fastq, Forward Sequence2.fastq, Forward Sequence3.fastq >> Merged_tissues_forward.fastq  
    cat Reverse Sequence1.fastq, Reverse Sequence2.fastq, Reverse Sequence3.fastq >> Merged_tissues_reverse.fastq  

#### 결과
**Merged_tissues_forward.fastq**  
**Merged_tissues_reverse.fastq**  

# **3) *De novo* assembly**
  
Trinity를 활용하여 assembly 진행
완료되면 여러 파일들이 생성되는데 이때  Trinity.fasta 라는 파일이 assembly 된 transcriptome  
이 과정 이후에 statistics 를 구하고 싶으면 Trinity tool 의 util directory 내의 TrinityStats.pl을 이용  


File: **Merged_tissues_forward.fastq, Merged_tissues_reverse.fastq**  
Tool: **Trinity**  

#### 커맨드

    Trinity --seqType fq --left Merged_tissues_1.fastq --right Merged_tissues_2.fastq --output trinity_out --max_memory 100G --CPU 8

#### 옵션
--seqType: reads format 을 지정 (fq: fastq, fa: fasta)  
--left, right: foward, reverse reads 를 지정  
--output: 생성될 output directory (trinity 단어가 들어가야함)  
--max_memory:  assembly 과정중 할당할 최대 memory  
--CPU:  사용할 cpu 수  

#### 결과
**Trinity.fasta**  
**TrinityStats.pl**  


# **4) Gene prediction**

## 4.1) 100 amino acids 이상의 ORF 포함하는 transcripts 추출
  
Assembled transcripts 의 protein coding genes 을 찾기 위해 TransDecoder 를 이용  

File: **Trinity.fasta**  
Tool: **TransDecoder**  

#### 커맨드

    TransDecoder.LongOrfs -t Trinity.fasta

#### 옵션
-t: assembly 를 통해 생성된 Trinity.fasta  

#### 결과
**longest_orfs.pep** (이후 분석에 사용)  

**longest_orfs.gff3**  

**longest_orfs.cds** 등 생성  

## 4.2) Homology-based search

기존에 잘알려진 protein database (Swissprot DB) 에 homology search 를 하여 gene prediction 과정에 활용  
  

### 4.2.1) swissprot.fasta 파일을 download 후 blast db 형식으로 제작

File: **sprot.fasta**  
Tool: **Blast**  

#### 커맨드

    ncbi-blast-2.x.x+/bin/makeblastdb -in swissprot.fasta -dbtype prot  

#### 옵션
-in:  db 형식으로 만들고자 하는 fasta 파일  
-dbtype:  fasta 파일이 nucleotide 또는 protein sequence 인지에 따라 (nucleotide: nucl,  protein: prot)  

#### 결과

**swissprot.fasta**  

---

### 4.2.2) db 파일에 transcript assembly 를 homology search

File: **logest_orfs.pep, swissprot.fasta**  
Tool: **Blast**  

#### 커맨드

    ncbi-blast-2.x.x+/bin/blastp -query longest_orfs.pep -db swissprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6  

#### 옵션
-query: db에 search 할 sequence 파일  

-db: database  

-max_target_seqs: 해당 query가 matching 되는 sequence 중 몇개를 보여줄지에 대한 parameter  

-outfmt: output의 format에 대한 정보 (0~11 까지 있으며 6 은 tabular format)  

-evalue: blast evalue cutoff  

-num_threads: 사용할 cpu 수  

#### 결과
**blastp.outfmt6**  

---

## 4.3) Gene prediction

이제 blast 결과를 토대로 gene prediction 을 진행  
transcripts.transdecoder_dir 가 있는 directory 에서 진행  

File: **Trinity.fasta, blastp.outfmt6**  
Tool: **TransDecoder**  

#### 커맨드

    TransDecoder.Predict -t Trinity.fasta --retain_blastp_hits blastp.outfmt6 --cpu 8 --single_best_only

#### 옵션
-t:  Trinity assembly 파일  

--retain_blastp_hits: (4-2-2) 과정에서 얻은 blastp 결과 파일  

--cpu: 사용할 cpu 수  

--single_best_orf: gene prediction 과정 중 하나의 transcripts 에서 6개의 ORFs 를 탐색했을때 가장 best orf 출력  

#### 결과
**Trinity.fasta.transdecoder.pep** 등 생성  

---

## 4.4) Removing redundant transcripts
  
(4-2-1)~(4-2-3) 과정을 통해 gene prediction 은 완료되었으나, transcriptome 이다 보니 불가피하게 isoforms 이 존재  
Reference genome 의 역할을 하기 위해 representative protein coding genes 만 있는게 좋으므로, redundant 한 transcripts 를 제거하기 위해 CD-hit 사용  
.cdhit 파일은 sequence,  .cdhit.clstr 파일에는 cluster 가 어떻게 묶였는지에 대한 정보가 담긺  
얻어진 Trinity.fasta.transdecoder.pep.cdhit 파일을 Non-redundant protein coding sequences (NRCDS) 로 명명하고 진행  

File: **Trinity.fasta.transdecoder.pep**  
Tool: **CDHit**  

#### 커맨드

    cdhit -i Trinity.fasta.transdecoder.pep -o Trinity.fasta.transdecoder.pep.cdhit -c 0.99 -T 8

#### 옵션
-i:  TransDecoder 를 통해 geneprediction 이 된 fasta 파일  

-o: otuput   

-c: identity cutoff (0.99: 99% 비슷한 transcripts 의 cluster 중 가장 긴 transcript 선택)  

-T:  사용할 cpu 수  

#### 결과
**Trinity.fasta.transdecoder.pep.cdhit**  
**Trinity.fasta.transdecoder.pep.cdhit.clstr**  

---

# 5) Gene expression level quantification

## 5.1) NRCDS_Trinity.fasta 파일 생성
  
NRCDS 파일에 대응되는 nucleotide sequence 가 필요, NRCDS 파일의 sequence id 를 *de novo* assembly 파일인 Trinity.fasta 에 matching 하여 nucleotide sequence 를 얻는다.  
    
커맨드 줄 해석  
1) ">" 가 있는 줄에서 첫번째 필드만 가져옴  
2) .p 제거  
3) .p 뒤에 있던 숫자 제거  
4) 앞의 > 제거  
5) TrinityID와 neucleotide sequence mapping  
  

File: **Trinity.fasta.transdecoder.pep.cdhit, Trinity.fasta**  
Tool: **Linux command, Samtools faidx**  

#### 커맨드

    awk '/>/ {print $1}' Trinity.fasta.transdecoder.pep.cdhit > output.txt
    sed 's/.p/ /g' output.txt > output2.txt
    awk '{print $1}' output2.txt > output3.txt
    sed 's/>//g' output3.txt > TrinitypepID.txt
    xargs samtools faidx Trinity.fasta <TrinitypepID.txt> NRCDS_sequence.fasta

#### 옵션
xargs: 대량으로 진행  
file1.fasta <[ID.txt]> file2.fasta: 괄호 안의 ID파일과 일치하는 file 1의 시퀀스 파일을 파싱해서 file2로 저장  


#### 결과
**NRCDS_sequence.fasta**  

## 5.2) Read alignment & Abundance estimation
  
이제 NRCDS_Trinity.fasta 에 reads 를 mapping 할 차례이다. 이는 Trinity 의 util 내의 align_and_estimate_abundance.pl script를 사용  
여기서의 reads는 concatenation된 reads가 아닌, 각각의 sample 의 reads  
결과 파일에는 각 transcripts 의 read count, TPM, FPKM 정보가 담겨있음  
각 sample 별로 분석을 진행  

File: **NRCDS_Trinity.fasta, rawdata_1.fastq, rawdata_2.fastq**  
Tool: **Trinity util 내의 util/align_and_estimate_abundance.pl**  

#### 커맨드

    Trinity tool path/util/align_and_estimate_abundance.pl --transcripts NRCDS_Trinity.fasta --seqType fq --left rawdata_1.fastq --right  rawdata_2.fastq --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir rawdata_rsem_out --output_prefix  rawdata_rsem --thread_count 8  

#### 옵션
--transcripts:  mapping 하고자 하는 assembly 파일 (NRCDS_Trinity.fasta)

--seqType:  reads format 을 지정 (fq: fastq, fa: fasta)

--left, right:  foward, reverse reads 를 지정

--est_method:  abundance estimation 하는 method 설정 (RSEM을 사용한다.)

--aln_method:  read alignment 하는 method 설정 (bowtie2를 사용한다.)

--prep_reference:  assembly 파일을 indexing 하는 과정 (처음에 한 sample 돌릴때 했다면 다음 sample 부터는 빼도 됨)

--output_dir: 생성될 output directory

--output_prefix: 생성될 output 파일 앞에 붙는 이름

--thread_count:  사용할 cpu 수  

#### 결과
**.RSEM.isoforms.results**  

# 6) DEGs(Differentially expressed genes) analysis  
  
각 sample 간의 DEGs 를 확인하고자 진행한다.  

## 6.1) Gene expression matrix
  
각 sample 로 부터 얻어진 RSEM.isoforms.results 를 하나의 matrix 로 합친다.  
각 sample 의 RSEM.isoforms.results 를 나열  

File: **rawdata1.RSEM.isoforms.results, rawdata2.RSEM.isoforms.results**  
Tool: **Trinity 내의 util/abundance_estimates_to_matrix.pl**

#### 커맨드

    Trinity tool path/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix ABC  rawdata1.RSEM.isoforms.results rawdata2.RSEM.isoforms.results -gene_trans_map none --name_sample_by_basedir  

#### 옵션
--est_method:  abundance estimation 할때 썻던 method  

--out_prefix: 생성될 output 파일 앞에 붙는 이름  

#### 결과
**(outprefix).counts.matrix를 포함한 여러 파일**  

## 6.2) Differential expression analysis  
  
각 sample의 발현값을 토대로 DE 분석을 진행한다.  
여러 sample 인 경우에 모든 1:1 조합으로 비교하기 때문에 n combination 2 의 경우의 수가 나온다.  
그러므로 많은 파일들이 생성됨  

File: **(outprefix).counts.matrix**  
Tool: **Trinity 내의 Analysis/DifferentialExpression/run_DE_analysis.pl**

#### 커맨드

    Trinity tool path/Analysis/DifferentialExpression/run_DE_analysis.pl  --matrix  counts.matrix  --method  edgeR  --output  edgeR_dir  --dispersion  0.1   

#### 옵션
--matrix:  앞서 얻었던 couts.matrix 파일  

--method: DE 분석에 사용할 method (R package 인 edgeR)  

--output:  output directory 이름  

--dispersion:  같은 종의 sample 인 경우 0.1 을 사용 (edgeR manual 에 자세한 설명이 있으니 참고)  

#### 결과
**많은 파일**  

## 6.3) TMM normalization
  
각 sample 에서 얻은 gene expression level 을 normalization 하는 과정  

### 6.3.1) trasncript length 정보 준비
  
1,3,4 번의 field (column) 만 추려내는 것으로 id, length, effective_length 정보를 준비한다.

File: **rawdata1.RSEM.isoforms.results**  
Tool: **Linux command**

#### 커맨드

    cut -f 1,3,4 A.RSEM.isoforms.results > Trinity.trans_lengths.txt    

#### 결과
**Trinity.trans_lengths.txt**  

---

### 6.3.2) Normalization  
  
RNA-seq normalization 방법중 하나인 TMM normalization 방법을 사용  

File: **isoform.counts.matrix**  
Tool: **Trinity 의 Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl**

#### 커맨드

    Trinity tool path/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix counts.matrix --length Trinity.trans_lengths.txt  

#### 옵션
--matrix:  앞서 얻었던 isoform.counts.matrix 파일  

--length:  앞서 얻었던 Trinity.trans_lengths.txt 파일  

#### 결과
**.counts.matrix.TMM_normalized.FPKM**  

## 6.4) PCA 분석 및 Heatmap 확인

PCA분석과 Heatmap을 이용하여 데이터의 분포를 확인  
이는 보통 같은 tissue인 샘플끼리 모이는 경향이 있는데 데이터가 제대로 나왔는지 확인도 하고 발현량이 다른 샘플을 확인하기 위해 사용  
Python으로 PCA를 하면 문제가 있다고 해서 R을 이용해 PCA rotation 데이터를 가져와 Python에서 그래프만 그림 <- 파이썬에서 제대로 안 되는것이 맞는지 확인필요  


### 6.4.1) R PCA rotation

Python PCA와 R PCA에는 차이가 좀 있는데 R PCA의 rotation 값이 필요하기 때문에 R PCA의 rotation 값을 가져와 PCA를 그림  
  
File: **Annotated TMM_FPKM**  
Tool: **R Code**  
  
#### 커맨드

    FPKM <- read.csv(File,header=True)
    FPKM <- FPKM[,c(2:sample_num)]
    FPKMpca <- prcomp(FPKM)
    write.csv(FPKMpca,"FPKMpca.csv")

### 6.4.2) PCA ,Heatmap 그래프 그리기

#### 커맨드

    R PCA,correlation.ipynb 참조

#### 결과
**PCA,Heatmap graph**  
  

## if ) Heatmap이 같은 조건의 샘플끼리 묶이지 않음 or 같은 샘플임에도 Correlation 값의 편차가 심하다면
### if. 1) Bacteria Contamination Check in Annotation Level

샘플, 조건 별로 TMM_normalized FPKM Matrix 값을 이용해 Correlation Graph를 그려 Gene 별 발현량 차이를 확인한다.  

비슷한 조건, 샘플임에도 특정 Gene의 발현량의 차이가 유난히 많이 보임 -> Bacteria Contamination Check  
  
샘플을 Blast로 Protein Annotation 한 결과로 확인을 하였는데 샘플을 Annotation한 Protein과 Bacteria가 만드는 Protein에 대한 DB와 비교해 보았다.  

Model organism이 아니기 때문인지는 몰라도 겹치는 Protein이 낮게 나오는 결과를 보였다.  

#### 커맨드

    PCA,correlation.ipynb 참조
    Bacteria Contam Check.ipynb 참조
  

### if. 2) Bacteria Contamination Check in Expression Level
  
전체 Bacteria Gene을 다운받아 cat하여 Bowtie2를 이용해 Trimming 된 샘플 데이터를 이용하여 Bacteria Gene에 얼마나 Mapping 되는지 확인한다.(Mapping 많이 될수록 많은 Contam)  
  
#### 커맨드

    1. NCBI -> Taxonomy -> Taxonomy Browser -> Bacteria -> Bacteria -> Genome subtree link -> Bacteria Genome download
    2. 압축 해제 후 리눅스 쉘에서 cat Bacteria_Genome -> Merged_bact.fasta
    3. bowtie2-build --threads 70 Merged_bact.fasta BacT
    4. bowtie2 --threads 10 -x BacT -1 3B1r_1_paired.fastq -2 3B1r_2_paired.fastq -S 3B1.sam

  
### if. 3) Clean Reads *De novo* Assembly 

Bowtie2 Mapping 과정에서 나온 SAM file을 samtools를 이용해 BAM file로 바꾸어준다.  

BAMfile에서 unmapped 된 시퀀스만 추출하여 Fastq파일로 만들어준다.  

샘플 별 Fastq파일을 cat 해준 뒤 Trinity 실행하여 Contam된 시퀀스를 제거한 clean 시퀀스로 De novo assembly를 진행한다.  

#### 커맨드

    1. samtools view -bS 3B1.sam > 3B1.bam
    2. samtools view -bf 0x04 out.bam > unmap3B1.bam
    3. samtools bam2fq unmap3B1.bam > clean3B1.fastq
    4. cat clean3B1.fastq | grep '^@.*/1$' -A 3 --no-group-separator > clean3B1_r1.fastq
    5. cat clean3B1.fastq | grep '^@.*/2$' -A 3 --no-group-separator > clean3B1_r2.fastq
    6. cat All_cleanSamples_r1.fastq > clean_merge_r1.fastq
    7. cat All_cleanSamples_r2.fastq > clean_merge_r2.fastq
    8. Trinity --seqType fq --left cleanMerge_r1.fastq --right cleanMerge_r2.fastq --output trinity_clean_dendro --max_memory 500G --CPU 40  

### 4번으로 돌아가서 다시 분석하시오.  
  
  

# 7) Annotation
  
얻어진 유전자들에 대해 기능이 잘 알려진 database 에 homology-based search 를 통해 유사한 기능을 할것이라고 추측  

사용할 수 있는 database 는 근연종이 genome 이 있다면 해당 genome data  
 
그렇지 않다면 swissprot, KEGG, UniRef 등에 BLAST 를 통해 해당 유전자들의 기능을 유추  

## 7.1) Homolgy Search

File: **swissprot.fasta, Trinity.fasta.transdecoder.pep.cdhit**  
Tool: **Blast**

#### 커맨드

    ncbi-blast-2.x.x+/bin/blastp -query Trinity.fasta.transdecoder.pep.cdhit -db swissprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6  

#### 옵션
-query: db에 search 할 sequence 파일  

-db: database  

-max_target_seqs: 해당 query가 matching 되는 sequence 중 몇개를 보여줄지에 대한 parameter  

-outfmt: output의 format에 대한 정보 (0~11 까지 있으며 6 은 tabular format)  

-evalue: blast evalue cutoff  

-num_threads: 사용할 cpu 수  

#### 결과
**blastp.outfmt6** 위의 blastp.outfmt6와 다른 파일임

## 7.2) Homolgy Search 후 중복 Annotated Genes 제거
  
Homolgy Search 후 중복 Annotated 유전자 제거
가장 긴 Gene만 남김  
커맨드 설명
1. 

File: ** - **  
Tool: ** - **

#### 커맨드  

    sed 's/.p/ /g' dendro_blastp.outfmt6 > output.txt
    sed 's/|/ /g' output.txt > output2.txt
    cat output2.txt | awk '{$2=""; print $0}' > output3.txt
    awk '{print $1}' output3.txt > dendro_blastp_ID.txt
    xargs samtools faidx /home/KHJ/KIM/Denovo/dendropanax/3.transDecoder/1.NRCDS/NRCDS/NRCDS_sequence.fasta <dendro_blastp_ID.txt> NRCDS_blast.fasta
    output3, NRCDS_blast2.fasta scp로 파일 가져오기
    NRCDS_blast_length.ipynb 참조
    


## 7.3) DEGs Annotation concat
  
6의 과정에서 얻은 DEGs 정보와 Annotation 정보를 합쳐줌  

TMM normalized FPKM에서 Homology Search 결과인 blastp.outfmt6 파일의 정보만 뽑아내는 과정

File: **blastp.outfmt6, RSEM.isoform.counts.matrix.TMM_normalized.FPKM**  
Tool: **Python Code**

#### 커맨드

    Annotatied TMM_FPKM.ipynb 참조  

#### 결과
**Annotated TMM_FPKM**  
  

# 8) DEGs Analysis

## 8.1) Identifying DEGs
  
위에서 얻은  .counts.matrix.TMM_normalized.FPKM 파일을  edgeR_dir 로 이동  
  
결과에는 위의 기준 (-C, -P) 에 충족하는 DEGs 들이 모여있음  
  
log2.centered.dat가 붙은 파일들이 있는데 이는 절대적인 발현값이 아닌 상대적인 발현값으로 계산한 정보  

File: **isoform.counts.matrix.TMM_normalized.FPKM**  
Tool: **Trinity 의 Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl**

#### 커맨드

    Trinity tool path/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix .counts.matrix.TMM_normalized.FPKM -C 2 -P 0.001    

#### 옵션
--matrix: 앞서 얻었던 counts.matrix.TMM_normalized.FPKM 파일

-C: Log2 foldchange 로 2는 발현값이 4배 높거나 4배 낮은 기준 (DEGs cutoff)

-P: 통계값  corrected p-value (FDR) 로 0.001 은  0.001 미만의 기준 (DEGs cutoff)

#### 결과
**diffExpr 의 prefix 가 붙은 파일** 


## 8.2) DEGs David 
  
DEGs Annotated 파일에서 Protein ID를 David DB에 넣어 어떤 기능을 하는지 확인
찾으려고 하는 기능들이 있는지 확인  

File: **DEGs Annotated list**  
Tool: **David DB**   

#### 커맨드

     DAVID (https://david.ncifcrf.gov/) 접속
     Start analysis
     Upload Gene list -> Select identifier -> Gene list -> Submit list
     List -> Select species -> annotate
     Gene Ontology(BP,CC,MF), KEGG pathway 확인
     

#### 결과
**Gene function**  

<script src="https://giscus.app/client.js"
        data-repo="watchback/watchback.github.io"
        data-repo-id="R_kgDOK8k5WQ"
        data-category="Announcements"
        data-category-id="DIC_kwDOK8k5Wc4Cb7Lj"
        data-mapping="og:title"
        data-strict="0"
        data-reactions-enabled="1"
        data-emit-metadata="0"
        data-input-position="bottom"
        data-theme="preferred_color_scheme"
        data-lang="ko"
        crossorigin="anonymous"
        async>
</script>