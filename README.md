# prj_edna_refdb
Generation of a Reference Database for environmental DNA (eDNA) analysis at INBO.

The current workflow includes equal amounts of manual/expert curation as automated scripts. 

In this repository we aim to develop and maintain code only within the context of **12S eDNA amplicon** analysis **at INBO**.  

### Requirements
- INBO google account with access to Input Data in [PRJ_eDNA_Refdb](https://drive.google.com/drive/folders/0AElXLysx42mqUk9PVA)
- Unix environment or Docker installation

## Workflow overview
**Input Data**

Input Data is prepared and maintained by the scientist(s), including:
- FASTA files of relevant DNA sequences that include 12S gene
- List of expected species (gsheet)
- List of sequences to filter (gsheet) 
- List of expected conflicts/multihits and how to resolve (gsheet)

### Workflow loop
Generation of a new reference database (on change of Input Data) can be executed by anyone.

This is done mainly in an **interactive R-studio session**:
1. **Pre-process** Input Data in R (i.e. apply filters)
2. Run **in-silico PCR** (not in R)
3. **Evaluate** output in R and google-sheets => by scientist, **update Input Data** if needed
4. If Input Data updated, start from  **1.** untill no more (unaccounted) conflicts


## Run Local
#### Documentation
- google docs: [Genereren referentiedatabanken](https://docs.google.com/document/d/1KGxwikBzEHPiCVs2te2lHzkE8tL66tyVldPJFDs6hKc/) 

## Run on Unix Server (HPC)
#### Documentation
- Set-up: [PRJ_eDNA_Refdb_2023-2025](https://www.milieuinfo.be/confluence/display/INBOBIO/PRJ_eDNA_Refdb_2023%3A2025)
- Interactive workflow: [[HOW-TO] Make an INBO-curated 12S reference database](https://www.milieuinfo.be/confluence/display/INBOBIO/%5BHOW-TO%5D+Make+a+INBO-curated+12S+reference+database)
### Prepare environment
- Launch an interactive Rstudio session via the [OnDemand web portal](https://ondemand.hpc.kuleuven.be/)
  - Don't forget to configure your rclone to mount the [PRJ_eDNA_Refdb](https://drive.google.com/drive/folders/0AElXLysx42mqUk9PVA) on the INBO gdrive
  - and load the `rclone` module via the pre-run scriplet
- Prepare a working directory on the system
- Download the (latest) NCBI taxdump

### Interactive workflow
- Run all `.R` scripts in chronological order
- Variable names in UPPERCASE might need to be provided as user-input
  - They are always declared in the beginning of each script
- When needed switch from `PRIMER_NAME="riaz"` to `PRIMER_NAME="teleo"`
- After pre-processing (`02`) you need to submit a job to run the `obi ecopcr` (done manually for now)
- After the job finishes (<5min) continue in Rstudio (`04`)

### Evaluation of created Reference Database

- Evaluation tables are written to corresponding google-sheets in the appropriate folder in [PRJ_eDNA_Refdb](https://drive.google.com/drive/folders/0AElXLysx42mqUk9PVA)
- These tables are discussed with the scientist(s)
- Any problematic conflicts are resolved by updating the input files on the google-drive [PRJ_eDNA_Refdb](https://drive.google.com/drive/folders/0AElXLysx42mqUk9PVA).
Some examples:
  - Add/remove sequences in the fasta files
  - Change/choose priority of target species in (multiple) google-sheets

### Output

- OBITools3 object in DMS format. Can be used directly with `obi ecotag`. See https://forge.metabarcoding.org/obitools/obitools3
- Exports of some VIEWS in fasta format to use with other taxonomic annotation tools

### Citation

Citation: Boyer F., Mercier C., Bonin A., Taberlet P., Coissac E. (2014) OBITools: a Unix-inspired software package for DNA metabarcoding. Molecular Ecology Resources, submitted.
