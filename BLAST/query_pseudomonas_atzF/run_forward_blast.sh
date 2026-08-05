# Add a header to the output file
printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n' > query_pseudomonas_atzF/results/forward_blast_hits.tsv

# Run BLASTp of atzF sequence against the RefSeq genome database
blastp -query query_pseudomonas_atzF/pseudomonas_atzF_seq.fa -db dbs/amac_refseq_db  -outfmt 6 -max_target_seqs 10 >> query_pseudomonas_atzF/results/forward_blast_hits.tsv