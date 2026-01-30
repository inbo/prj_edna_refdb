
source("scripts/refdb_teleo_00_initialisation.R")
source("scripts/_functions_fasta.R")


files <- sort(list.files(fasta_inputs, pattern = ".fasta"))

#read the source files
df_inputs_all <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(fasta_inputs, file))
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
}

origseqs <- read_csv2("database/refdb_teleo_2023-08-10/te_controleren_seqs.csv")

corrections_taxon <- df_inputs_all %>% 
  select(genlab_id, taxid_orig = taxid, dna_sequence) %>% 
  inner_join(origseqs %>% select(genlab_id, taxid = gevonden_taxon), by = "genlab_id") %>% 
  mutate(genlab_id = paste0(substring(genlab_id, 1, 15), "_MOD"))



create_input_fasta(file = "corrected_genlab_records.fasta", data = corrections_taxon)


### Kijken of de primer te vinden is

mtch <- NULL
for (i in seq_len(nrow(df_inputs)))
  mtch[i] <- find_matching_seq(df_inputs$dna_sequence[i])

test <- create_random_seq(size = 100, preamble = 5, postamble = 5)
test <- toupper("caaaggcttggtcctgactttactatcaactctagctaaacttacacatgcaagtatccgcatccccgtgagaatgccctacagttccctgcccgggaacaaggagctggtatcaggcacacttcgactagcccatgacaccttgcttagccacacccccaagggaactcagcagtgatagacattaagccataagtgaaaacttgacttagtcaaagctaagagggccggtaaaactcgtgccagccaccgcggttatacgagaggcccaagttgatagacatcggcgtaaagcgtggttaagattaaagacaatactaaagccgaacaccttcagagctgttatacgcatccgaaggtaagaagttcaaccacgaaagtggctttatagcccctgaacccacgaaagctacgatacaaactgggattagataccccactatgcctagccataaacattggtagcacactacacccactacccgcctgggaactacgagcatcagcttgaaacccaaaggacttggcggtgctttagatccacctagaggagcctgttctagaaccgataacccccgttcaacctcacctttccttgtctctcccgcctatataccgccgtcgtcagcttaccctgtgaaggttaaatagtaagcaaaattggtacaacctaaaacgtcaggtcgaggtgtagcgtatgggaagggaagaaatgggctacatttcctattacaggaaatacgaatggtgtactgaaacgtacgcctgaaggaggatttagcagtaagcaggaaatagagcgtcccgctgaaattggccctgaagcgcgcacacaccgcccgtcactctccccaagcctaccaactaaaataattaaaaccctataatcgcgaaggggaggcaagtcgtaacatggtaagtgtaccggaaggtgcacttggaaaaat")
find_matching_seq(test)

niet_ampli <- read_csv2("C:\\_GIT_PROJECTS\\prj_edna_refdb\\database\\refdb_teleo_2023-08-10\\output\\sequenties_niet_geamplificeerd.csv")

mtch <- NULL
for (i in seq_len(nrow(niet_ampli)))
  mtch[i] <- find_matching_seq(niet_ampli$dna_sequence[i])
table(mtch)

niet_ampli$mtch <- mtch

geweigerd_ondanks_match <- niet_ampli %>% filter(substring(mtch, 1,3) == "1_1") 

write_excel_csv2(geweigerd_ondanks_match, file = file.path(output_path, "geweigerd_ondanks_perfecte_match.csv"))








