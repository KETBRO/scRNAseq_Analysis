## For cellranger-9.0.1
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1746058118&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=QfF~h4-k6hzlLzga6WdIdoG7mXtBRQJKT4BVISEnsy5IODtzDDU8uUjkzE3U-O-LDQMRsRQFG9Lq-Z8uuHaBS6Z-VBrJT42C4cqq6dWSPyphhrR87raWNxxFz~pBT6rW7LJp6JyayJmlR0N4nNYJglg6ACcw7r-v5kForTnaa8vhvZrmMYr8B2pxXCsQBXgSt26XjjZknuGheiHB6utkOFR4zeVttVaejxnzCtrs8fk2iZWqkmCuCgqtWjX3IDfW5TfpffAIzQCKypWFY0yDyoTfpZWOyWqTpR9JVWFSgnHKFLrBvrTGFWhlrP8gt3PCpmA~JxPZITaP0A9-adeIBA__"

tar -xvzf cellranger-9.0.1.tar.gz
sudo mv cellranger-9.0.1 /opt/
sudo ln -s /opt/cellranger-9.0.1/cellranger /usr/local/bin/cellranger

#check cell ranger version
cellranger --version

#Extract FASTQ files
mkdir 3p_Citrate_CPT_fastqs
tar -xvf 3p_Citrate_CPT_fastqs.tar -C 3p_Citrate_CPT_fastqs

# Download or generate Cell Ranger reference
# Human GRCh38 reference from 10x
wget -O refdata-gex-GRCh38-2020-A.tar.gz "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
tar -xvzf refdata-gex-GRCh38-2020-A.tar.gz

#Run Cell Ranger Count
cellranger count \
  --id=citrate_sample \
  --transcriptome=refdata-gex-GRCh38-2020-A \
  --fastqs=3p_Citrate_CPT_fastqs \
  --sample=sample \
  --localcores=12 \
  --localmem=16

