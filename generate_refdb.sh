#!/bin/bash
cd /app
source /app/obi3-env/bin/activate
echo 'aantal inputs: ' > logfile.txt
grep -E  -i '>' input.fasta | wc -l >> logfile.txt
obi import input.fasta refdb/input
obi export refdb/input -o kept_input.fasta
echo 'aantal behouden inputs: ' >> logfile.txt
grep -E  -i '>'  kept_input.fasta | wc -l >> logfile.txt
obi import --taxdump taxdump.tar.gz refdb/taxonomy/dump
obi ecopcr -e 4 -l 50 -L 160 -F ACTGGGATTAGATACCCC -R TAGAACAGGCTCCTCTAG --taxonomy refdb/taxonomy/dump refdb/input refdb/ecopcr
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy refdb/taxonomy/dump refdb/ecopcr refdb/ecopcr_clean
obi export refdb/ecopcr_clean -o amplified_clean.fasta
echo 'aantal geamplificeerde inputs: ' >> logfile.txt
grep -E  -i '>'  amplified_clean.fasta | wc -l >> logfile.txt
obi uniq --taxonomy refdb/taxonomy/dump refdb/ecopcr_clean refdb/ecopcr_uniq
obi grep --require-rank=family --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq  refdb/ecopcr_uniq_clean
obi export refdb/ecopcr_uniq_clean -o amplified_clean_uniq.fasta
obi build_ref_db -t 0.97 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq_clean refdb/ecopcr_final_0.97
obi build_ref_db -t 0.99 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq_clean refdb/ecopcr_final_0.99
obi build_ref_db -t 1 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq_clean refdb/ecopcr_final_1
obi export  --fasta-output refdb/ecopcr_final_0.97 -o final_db_0.97.fasta
obi export  --fasta-output refdb/ecopcr_final_0.99 -o final_db_0.99.fasta
obi export  --fasta-output refdb/ecopcr_final_1 -o final_db_1.fasta
echo 'aantal unieke inputs na merging: ' >> logfile.txt
grep -E  -i '>'  final_db_0.97.fasta | wc -l >> logfile.txt

