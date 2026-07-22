# Add a header to the output file
printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n' > query_opa_family/results/forward_blast_hits.tsv

# Run BLAST search against the Amac database
blastp -query query_opa_family/opa_family_seqs.fa -db dbs/amac_db -max_target_seqs 10 -outfmt 6 >> query_opa_family/results/forward_blast_hits.tsv