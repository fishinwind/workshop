# from: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30227

extract() {
    f=$1
    echo $f
    #./bin/fastq-dump -M 40 --read-filter pass --qual-filter --gzip $f
    fastq-dump --origfmt --qual-filter --gzip $f
}

# GM12878_DNAse_REP1
wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX017%2FSRX017008/SRR036667/SRR036667.sra"
extract ./SRR036667.sra
mv SRR036667.fastq.gz GM12878_DNAse_REP1.fastq.gz

# GM12878_DNAse_REP2
wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX017%2FSRX017009/SRR036668/SRR036668.sra"
extract ./SRR036668.sra
mv SRR036668.fastq.gz GM12878_DNAse_REP2.fastq.gz

# GM12878_DNAse_REP3
wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX017%2FSRX017010/SRR036669/SRR036669.sra"
extract ./SRR036669.sra
mv SRR036669.fastq.gz GM12878_DNAse_REP3.fastq.gz


# GM12878_input DNA
wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX080%2FSRX080134/SRR298996/SRR298996.sra"
extract ./SRR298996.sra
mv SRR298996.fastq.gz GM12878_input.fastq.gz

