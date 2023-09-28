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

obi ecopcr -e 5 -l 20 -L 100 -F ACACCGCCCGTCACTCT -R CTTCCGGTACACTTACCRTG --taxonomy refdb/taxonomy/dump refdb/input refdb/ecopcr
obi export refdb/ecopcr -o amplified.fasta
echo 'aantal geamplificeerde inputs: ' >> logfile.txt
grep -E  -i '>'  amplified.fasta | wc -l >> logfile.txt

#cleaning:
#obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy refdb/taxonomy/dump refdb/ecopcr refdb/ecopcr_clean
#obi export refdb/ecopcr_clean -o amplified_clean.fasta
#echo 'aantal gegrepte geamplificeerde inputs: ' >> logfile.txt
#grep -E  -i '>'  amplified_clean.fasta | wc -l >> logfile.txt

#obi uniq --taxonomy refdb/taxonomy/dump refdb/ecopcr_clean refdb/ecopcr_uniq #indien cleaning
obi uniq --taxonomy refdb/taxonomy/dump refdb/ecopcr refdb/ecopcr_uniq #geen cleaning

obi export refdb/ecopcr_uniq -o amplified_uniq.fasta
echo 'aantal uniek geamplificeerde inputs: ' >> logfile.txt
grep -E  -i '>'  amplified_uniq.fasta | wc -l >> logfile.txt

#cleaning
#obi grep --require-rank=family --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq  refdb/ecopcr_uniq_clean 
#obi export refdb/ecopcr_uniq_clean -o amplified_clean_uniq.fasta 
#obi export refdb/ecopcr_uniq_clean -o amplified_uniq.fasta
#echo 'aantal uniek met grep gecleande inputs: ' >> logfile.txt
#grep -E  -i '>'  amplified_clean_uniq.fasta | wc -l >> logfile.txt

#obi build_ref_db -t 0.97 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq_clean refdb/ecopcr_final_0.97 #indien cleaning
#obi build_ref_db -t 0.995 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq_clean refdb/ecopcr_final_0.99 #indien cleaning
obi build_ref_db -t 0.97 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq refdb/ecopcr_final_0.97 #geen cleaning
obi build_ref_db -t 0.995 --taxonomy refdb/taxonomy/dump refdb/ecopcr_uniq refdb/ecopcr_final_0.99 #geen cleaning

obi export  --fasta-output refdb/ecopcr_final_0.97 -o final_db_0.97.fasta
obi export  --fasta-output refdb/ecopcr_final_0.99 -o final_db_0.99.fasta

echo 'aantal unieke inputs na merging: ' >> logfile.txt
grep -E  -i '>'  final_db_0.97.fasta | wc -l >> logfile.txt
grep -E  -i '>'  final_db_0.99.fasta | wc -l >> logfile.txt

