#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch 
#SBATCH --account=lp_inbo_pilot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:15:00
#SBATCH --out=obi-refdb-12S-%j.out
#SBATCH --job-name=obi_refdb_12S_iteration

#### INPUT:
# LAUNCH DIR: database/
# - IN_FASTA: input sequences cleaned with Rscripts prj_eDNA_REFDB (FASTA)
# - TAX_DUMP_GZ: zipped ncbi taxdump
# - PRIMER_NAME: type primer: riaz or teleo
#### EXAMPLE job submission:
# sbatch .../refdb.slurm -i input.fasta -t taxdump.tar.gz -p teleo

################################
# Commandline Argument Parsing #
################################

# default
PRIMER_NAME="riaz"

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      IN_FASTA=$2
      shift 2;;
    -t|--taxdmp)
      TAX_DUMP_GZ=$2
      shift 2;;
    -p|--primer)
      PRIMER_NAME=$2
      shift 2;;
    -*|--*)
      echo unknown argument $1
      exit 1;;
  esac
done

# Set primer specific parameters
#Riaz or Teleo analysis?
if [ $PRIMER_NAME == "riaz" ]; then
  echo "RIAZ-ecopcr"
  PRIMER_MINL=50
  PRIMER_MAXL=160
  PRIMER_E=4
  PRIMER_FWD=ACTGGGATTAGATACCCC
  PRIMER_REV=TAGAACAGGCTCCTCTAG

elif [ $PRIMER_NAME == "teleo" ]; then
  echo "TELEO-ecopcr"
  PRIMER_MINL=20
  PRIMER_MAXL=200
  PRIMER_E=5
  PRIMER_FWD=ACACCGCCCGTCACTCT
  PRIMER_REV=CTTCCGGTACACTTACCRTG

else
	exit 1
fi

###############
# PREPARATION #
###############

module purge
module load Python/3.9.5-GCCcore-10.3.0
source obitools3-venv

# INIT Log
MY_LOG="$SLURM_JOB_ID"_logfile.txt
echo 'aantal inputs: ' > $MY_LOG
grep -c "^>" $IN_FASTA >> $MY_LOG

######################
### OBITOOLS REFDB ###
######################

OBI_REF_DB=$(date -I)-obitools3-refdb-"$PRIMER_NAME"

rm -r "$OBI_REF_DB";mkdir "$OBI_REF_DB"

OBI_REF_DB="$OBI_REF_DB/$OBI_REF_DB"

obi import $IN_FASTA $OBI_REF_DB/input

obi export $OBI_REF_DB/input -o "$OBI_REF_DB"-input.fasta

echo 'aantal behouden inputs: ' >> $MY_LOG
grep -c "^>"  "$OBI_REF_DB"-input.fasta >> $MY_LOG

obi import --taxdump $TAX_DUMP_GZ $OBI_REF_DB/taxonomy/taxdump

obi ecopcr -e $PRIMER_E -l $PRIMER_MINL -L $PRIMER_MAXL \
  -F $PRIMER_FWD -R $PRIMER_REV \
  --taxonomy $OBI_REF_DB/taxonomy/taxdump $OBI_REF_DB/input $OBI_REF_DB/ecopcr

obi export $OBI_REF_DB/ecopcr -o "$OBI_REF_DB"-ecopcr.fasta

echo 'aantal geamplificeerde inputs: ' >> $MY_LOG
grep -c "^>"  "$OBI_REF_DB"-ecopcr.fasta >> $MY_LOG

#cleaning:
#obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy refdb/taxonomy/taxdump refdb/ecopcr refdb/ecopcr_clean
#obi export refdb/ecopcr_clean -o amplified_clean.fasta
#echo 'aantal gegrepte geamplificeerde inputs: ' >> logfile.txt
#grep -E  -i '>'  amplified_clean.fasta | wc -l >> logfile.txt

#obi uniq --taxonomy refdb/taxonomy/taxdump refdb/ecopcr_clean refdb/ecopcr_uniq #indien cleaning
obi uniq --taxonomy $OBI_REF_DB/taxonomy/taxdump $OBI_REF_DB/ecopcr $OBI_REF_DB/ecopcr_uniq #geen cleaning

obi export $OBI_REF_DB/ecopcr_uniq -o $OBI_REF_DB-ecopcr_uniq.fasta

echo 'aantal uniek geamplificeerde inputs: ' >> $MY_LOG
grep -c "^>" $OBI_REF_DB-ecopcr_uniq.fasta >> $MY_LOG

#cleaning
#obi grep --require-rank=family --taxonomy refdb/taxonomy/taxdump refdb/ecopcr_uniq  refdb/ecopcr_uniq_clean 
#obi export refdb/ecopcr_uniq_clean -o amplified_clean_uniq.fasta 
#obi export refdb/ecopcr_uniq_clean -o amplified_uniq.fasta
#echo 'aantal uniek met grep gecleande inputs: ' >> logfile.txt
#grep -E  -i '>'  amplified_clean_uniq.fasta | wc -l >> logfile.txt

#obi build_ref_db -t 0.97 --taxonomy refdb/taxonomy/taxdump refdb/ecopcr_uniq_clean refdb/ecopcr_final_0.97 #indien cleaning
#obi build_ref_db -t 0.995 --taxonomy refdb/taxonomy/taxdump refdb/ecopcr_uniq_clean refdb/ecopcr_final_0.99 #indien cleaning
obi build_ref_db -t 0.97 --taxonomy $OBI_REF_DB/taxonomy/taxdump $OBI_REF_DB/ecopcr_uniq $OBI_REF_DB/ecopcr_final_0.97 #geen cleaning
obi build_ref_db -t 0.995 --taxonomy $OBI_REF_DB/taxonomy/taxdump $OBI_REF_DB/ecopcr_uniq $OBI_REF_DB/ecopcr_final_0.995 #geen cleaning

obi export  --fasta-output $OBI_REF_DB/ecopcr_final_0.97 -o $OBI_REF_DB-ecopcr_final_0.97.fasta
obi export  --fasta-output $OBI_REF_DB/ecopcr_final_0.995 -o $OBI_REF_DB-ecopcr_final_0.995.fasta

echo 'aantal unieke inputs na merging: ' >> $MY_LOG
grep -c "^>" $OBI_REF_DB-ecopcr_final_0.97.fasta >> $MY_LOG
grep -c "^>" $OBI_REF_DB-ecopcr_final_0.995.fasta >> $MY_LOG

# Final obi refdb object
obi ls  $OBI_REF_DB
