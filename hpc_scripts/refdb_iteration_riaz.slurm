#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch 
#SBATCH --account=lp_inbo_pilot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:15:00
#SBATCH --out=obi-refdb-riaz-%j.out
#SBATCH --job-name=obi_refdb_riaz_iteration

#### INPUT:
# LAUNCH DIR: database/
# - IN_FASTA: input sequences cleaned with Rscripts prj_eDNA_REFDB (FASTA)
# - TAX_DUMP_GZ: zipped ncbi taxdump
# - PRIMER_NAME: type primer: riaz or teleo
#### EXAMPLE job submission:
# sbatch .../refdb.slurm -i input.fasta -t taxdump.tar.gz -p teleo

