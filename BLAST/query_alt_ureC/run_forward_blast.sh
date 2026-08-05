# Add a header to the output file
printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n' > query_alt_ureC/results/forward_blast_hits.tsv

# Run BLASTp of PG1 dut sequence against the RefSeq genome database
blastp -query query_alt_ureC/alt_ureC_seqs.fa -db dbs/amac_refseq_db  -outfmt 6 -max_target_seqs 10 >> query_alt_ureC/results/forward_blast_hits.tsv