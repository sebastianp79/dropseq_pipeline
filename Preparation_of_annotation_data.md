##Preparation of anotation files for the Drop-seq Pipelines



### 1. Obtain reference genomes



hg38
```bash
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz .
#combine into one single drop alternative haplotypes (_alt)
tar xvzf hg38.chromFa.tar.gz
cd chroms
rm *_alt.fa
cat *.fa > ../hg38_UCSC/hg38_ucsc.fa
rm -r chroms
```



mm10

```bash
#chromFa.tar.gz - The assembly sequence in one file per chromosome.
#    Repeats from RepeatMasker and Tandem Repeats Finder (with period
#    of 12 or less) are shown in lower case; non-repeating sequence is
#    shown in upper case.

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz .
tar xvzf chromFa.tar.gz
cat *.fa > mm10_UCSC/mm10_ucsc.fa
rm *.fa

```

danRer10 (zebrafish)

```bash
#Danio ri.. zebrafish

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz .

gunzip danRer10.fa.gz

#the model system expresses additional markers that have to be included in the reference genome:

EGFP


mCherry



mouse c-myc

```

panTro5

```bash




```



### 2. Obtain transcript annotation

hg38
```bash
Release 27 (GRCh38.p10)
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
```


mm10
```bash
Release M15 (GRCm38.p5)
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.annotation.gtf.gz

```

```bash
# downloaded danRer10 from UCSC
RefGene_danRer10-2017-08-28.gtf

#in addition add the annotations for the reporter genes EGFP and mCherry, and the overexpressed oncogene c-myc




```
