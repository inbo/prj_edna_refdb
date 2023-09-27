## STAP 1: ZORG DAT DOCKER DESKTOP RUNT en OBITOOLS geladen is
docker_container <- "obitools"

# ECOPCR

## instelvariabelen

e <- 4
l <- 50
L <- 160
p1 <- select_primer_tags("RIAZ")["oligo1"]
p2 <- select_primer_tags("RIAZ")['oligo2']
db_namepath <- paste0(db_name, "/", db_name) 

## maak structuur aan en kopieer de input.fasta naar docker
system2("docker", args = c("exec", docker_container, "mkdir", db_name))
system2("docker", args = c("cp", file.path("database", db_name, fasta_name), 
                           paste0(docker_container,":/",db_name)))

## converteer input+taxonomie naar voor obitools bruikbaar formaat
system2("docker", 
        args = c("exec", docker_container, "obiconvert", "--fasta", 
                 paste0("--ecopcrdb-output=",db_namepath),
                 "-t", taxdump_name, 
                 paste0(db_name, "/", fasta_name)))

## voer de in silico PCR uit
ecopcrstring <- paste("ecoPCR", "-d", db_namepath, 
                      "-e", e, "-l", l, "-L", L, p1, p2, 
                      ">", paste0(db_name, "/", "amplified.ecopcr"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', 
                           paste0('"', ecopcrstring, '"')))

## clean:(hou enkel die over die een parent rank hebben die soort en familie zijn)

grepstring <- paste("obigrep", "-d", db_namepath, 
                    "--require-rank=species",
                    "--require-rank=family",
                    paste0(db_name, "/", "amplified.ecopcr"),
                    ">", 
                    paste0(db_name, "/", "amplified_clean.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', 
                           paste0('"', grepstring, '"')))

## Keep unique sequences (merge identical amplicon sequences)

uniqstring <- paste("obiuniq", "-d", db_namepath, 
                    paste0(db_name, "/", "amplified_clean.fasta"),
                    ">", 
                    paste0(db_name, "/", "amplified_clean_uniq.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', 
                           paste0('"', uniqstring, '"')))

## Copy back to windows
system2("docker", args = c("cp", paste0(docker_container, ':', db_name),
                           file.path("database")))
