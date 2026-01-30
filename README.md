# prj_edna_refdb
Creation of Reference Database for E-DNA analysis at INBO

## Documentation
Internal documentation:
- google docs: [Genereren referentiedatabanken](https://docs.google.com/document/d/1KGxwikBzEHPiCVs2te2lHzkE8tL66tyVldPJFDs6hKc/)
- confluence: [PRJ_eDNA_Refdb_2023-2025](https://www.milieuinfo.be/confluence/display/INBOBIO/PRJ_eDNA_Refdb_2023%3A2025)

## Workflow (HPC)
- See also: [[HOW-TO] Make a INBO-curated 12S reference database](https://www.milieuinfo.be/confluence/display/INBOBIO/%5BHOW-TO%5D+Make+a+INBO-curated+12S+reference+database)
### Prepare environment
- Launch an interactive Rstudio session via the [OnDemand web portal](https://ondemand.hpc.kuleuven.be/)
  - Dont forget to first configure your rclone to mount the input-data-folder on the INBO gdrive ([PRJ_eDNA_Refdb](https://drive.google.com/drive/folders/0AElXLysx42mqUk9PVA))
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
