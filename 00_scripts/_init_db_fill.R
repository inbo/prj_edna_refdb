library(tidyverse)
library(DBI)
library(RSQLite)
source("scripts/_functions_fasta.R")


refdb_name <- "ReferentieDB.sqlite"
mydb <- RSQLite::dbConnect(RSQLite::SQLite(), refdb_name)

# INPUTS VAN REFDB 2022-05-17 (let op, ecopcr kan max 20 characters voor ID)

df_inputs <- parse_refdb_fasta("database/2022-05-17/input.fasta") %>% 
  transmute(ENTRY_ID = substring(genlab_id, 1, 20),
            SOURCE_DB = ifelse(substring(ENTRY_ID,1,1) == 1, "INBO", "GENBANK"),
            SOURCE_FILE = "",
            TAXID = as.numeric(taxid), 
            PRIMER = "RIAZ",
            LEN = nchar(dna_sequence),
            GENE = "12S",
            DNA_SEQ = dna_sequence,
            INPUT_HASH = AMPLICON_HASH,
            ACTIVE = 1)

DBI::dbAppendTable(mydb, name = "INPUT", df_inputs)

#error sequence (maar komt al niet meer voor in de inputs)
DBI::dbSendQuery(mydb, "update INPUT set ACTIVE = 0, REMARK = 'ERROR in TAXID' where ENTRY_ID = 'JN627425'")

# SOORTENLIJST VAN REFDB 2022-05-17

df_species_old <- read_csv2("database/2022-05-17/soortenevaluatie.csv") %>% 
  transmute(GRP, TAXID, RANK = INBO_RANK, 
            SCI_NAME, DUTCH_NAME, ENG_NAME, ACTIVE = 1,
            PRIORITY) %>% 
  distinct() %>% 
  filter(!(TAXID == 70285 & is.na(DUTCH_NAME)))
DBI::dbAppendTable(mydb, name = "SPECIES_LIST", 
                   df_species_old %>% select(-PRIORITY))

# SOORTENPRIORITEIT

project <- "RIAZ_GENERIC"
df_priority <- df_species_old %>% 
  transmute(TAXID, PRIORITY, PROJECT = project, ACTIVE = 1)
DBI::dbAppendTable(mydb, name = "SPECIES_PRIORITY", df_priority)

# ECOPCR

df_refdb_version <- data.frame(NAME = "RIAZ_2022-05-17",
                               DATE = "2022-05-17",
                               PRIMER = "RIAZ",
                               DESCRIPTION = "Referentiedb INBO Riaz 17 mei 2022")
DBI::dbAppendTable(mydb, name = "REFDB_VERSION", df_refdb_version)

df_ecopcr_settings <- readLines("database/2022-05-17/amplified.ecopcr", n = 13)
settings <- c("program", "version", 
              "oligo1", "oligo2", 
              "oligo2rc", "oligo1rc",
              "oligo1_Tmelt", "oligo2_Tmelt",
              "max_error_count", "dbname", "ampliconlength", 
              "taxdump_date")
settingvalues <- c("ecopcr-v2", "1.0.1", 
                   "ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG",
                   "CTAGAGGAGCCTGTTCTA", "GGGGTATCTAATCCCAGT",
                   "51.35", "51.12",
                   "4", "Riaz_2022_05_17", "[50,160] bp",
                   "2022-01-07"
                   )
df_ecopcr_settings <- data.frame(REFDB_VERSION_ID = 1, 
                                 SETTING = settings,
                                 VALUE = settingvalues)
DBI::dbAppendTable(mydb, name = "ECOPCR_SETTINGS", df_ecopcr_settings)

df_ecopcr <- read_delim("database/2022-05-17/amplified.ecopcr", 
                        delim = " | ", skip = 13, 
                        col_names = FALSE, trim_ws = TRUE) %>% 
  transmute(REFDB_VERSION_ID = 1, ENTRY_ID = X1, STRAND = X13, 
            SEQ_LEN_INPUT = X2, SEQ_LEN_AMPLICON = X20,
            AMPLICON = X21, AMPLICON_HASH = sha1hash(AMPLICON),
            RANK = X4, TAXID = X3, 
            SPECIES = X5, SPECIES_NAME = X6,
            GENUS = X7, GENUS_NAME = X8, 
            FAMILY = X9, FAMILY_NAME = X10, 
            FORWARD_MATCH = X14, FORWARD_ERROR = X15, FORWARD_TM = X16,
            REVERSE_MATCH = X17, REVERSE_ERROR = X18, REVERSE_TM = X19
            )
#nu nog IS_MERGED, OBI_COUNTm OBI_RANK, OBI_TAXID, MERGED_OVERVIEW uit processing
df_processed_ecopcr <- read_csv2("database/2022-05-17/ECOPCR_OUTPUT_SQLITE_2022-05-17.csv")

df_ecopcr <- df_ecopcr %>% 
  left_join(df_processed_ecopcr %>% 
              select(ENTRY_ID, IS_MERGED, OBI_COUNT, 
                     OBI_RANK, OBI_TAXID, MERGED_OVERVIEW))
#31 sequenties verloren gegaan bij de join (waar is_merged NA is)
table(is.na(df_ecopcr$IS_MERGED))
test <- df_ecopcr %>% filter(is.na(IS_MERGED))
write_excel_csv2(test, file = "technisch_teloorgegaan.csv")
#dit zijn ecopcr-fouten, omdat de taxonomie niet correct is toegewezen (taxid was genus, of genus/familie werden niet correct aan de soort toegewezen en zijn missing)

DBI::dbAppendTable(mydb, name = "ECOPCR_OUTPUT", df_ecopcr)


# MULTIHIT + CURATED + allowed higher level

## MULTIHIT
df_multihit <- read_csv2(file = "database/2022-05-17/multihit_en_blocklist.csv") %>% 
  transmute(TAXID, PREF_TAXID, PROJECT = "RIAZ_GENERIC", PRIMER = "RIAZ",
            EVALUATION = NA, BLOCKED_ON, ACTIVE = 1)

DBI::dbAppendTable(mydb, name = "MULTIHIT_TAXA", df_multihit)

## CURATED
curated_taxa <- df_multihit %>% 
  filter (TAXID == PREF_TAXID)  %>% 
  select(TAXID)
curated_seqs <- df_inputs %>% select(ENTRY_ID, TAXID) %>% 
  inner_join(curated_taxa) %>% 
  mutate(PROJECT = "RIAZ_GENERIC", ACTIVE = 1)

DBI::dbAppendTable(mydb, name = "CURATED_SEQS", curated_seqs)

## Allow reporting on level higher than species

df_allowed_merge <- 
  data.frame(TAXID = c(8033, 8041, 27772, 85421, 437164),
             RANK = "genus",
             PRIMER = "RIAZ",
             PROJECT = "RIAZ_GENERIC",
             SCI_NAME = c("Salvelinus", "Hucho", "Coregonus",
                          "Knipowitschia", "Achondrostoma"),
             BLAST_NAME = "bony fishes",
             ACTIVE = 1, 
             REMARK = NA)
DBI::dbAppendTable(mydb, name = "ALLOWED_MERGES", df_allowed_merge)



# VERSIECONTROLE
  # -> vul de ids in van alle tabellen die gebruikt zijn
 #input, species, priority, ecopcr, multihit, curated

inputids <- dbGetQuery(mydb, "select ID from INPUT") 
speciesids <- dbGetQuery(mydb, "select ID from SPECIES_LIST") 
priorityids <- dbGetQuery(mydb, "select ID from SPECIES_PRIORITY") 
ecopcrids <- dbGetQuery(mydb, "select ID from ECOPCR_OUTPUT") 
multihitids <- dbGetQuery(mydb, "select ID from MULTIHIT_TAXA") 
curatedids <- dbGetQuery(mydb, "select ID from CURATED_SEQS") 
alloweds <- dbGetQuery(mydb, "select ID from ALLOWED_MERGES") 

df_version_control <- 
  tibble(TABLE_NAME = "INPUT", 
         TABLE_ID = inputids %>% pull(ID)) %>% 
  bind_rows(tibble(TABLE_NAME = "SPECIES_LIST",
                   TABLE_ID = speciesids %>% pull(ID))) %>% 
  bind_rows(tibble(TABLE_NAME = "SPECIES_PRIORITY",
                   TABLE_ID = priorityids %>% pull(ID))) %>% 
  bind_rows(tibble(TABLE_NAME = "ECOPCR_OUTPUT",
                   TABLE_ID = ecopcrids %>% pull(ID))) %>% 
  bind_rows(tibble(TABLE_NAME = "MULTIHIT_TAXA",
                   TABLE_ID = multihitids %>% pull(ID))) %>% 
  bind_rows(tibble(TABLE_NAME = "CURATED_SEQS",
                   TABLE_ID = curatedids %>% pull(ID))) %>% 
  bind_rows(tibble(TABLE_NAME = "ALLOWED_MERGES",
                   TABLE_ID = alloweds %>% pull(ID))) %>% 
  mutate(REFDB_VERSION_ID = 1)

dbAppendTable(mydb, "REFDB_VERSION_CONTROL", df_version_control)


