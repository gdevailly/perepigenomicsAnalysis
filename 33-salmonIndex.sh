# Salmon
cd ~/work/projects/cascade/data/annotations/
salmon index -v
# salmon 0.9.1

salmon index --gencode -t gencode.v29.transcripts.fa.gz -i gencode.v29.transcripts.salmon.index

salmon quant
# --seqBias
# --gcBias
# -p 8
# -g gencode22.gtf
# --biasSpeedSamp -5

# Cool, it's working :)  v
salmon quant -i annotation/gencode.v24.transcripts.salmon.index -l A \
    -r rnaseq_reads/SRR020287.fastq.gz \
    -p 10 -g annotation/gencode.v24.annotation.gff3 \
    --seqBias --gcBias --biasSpeedSamp 5 \
    -o salmonQuant/


