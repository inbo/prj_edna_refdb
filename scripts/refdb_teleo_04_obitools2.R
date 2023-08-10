## STAP 1: ZORG DAT DOCKER DESKTOP RUNT en OBITOOLS geladen is
#run de container obitools_pv 

docker_container <- "obitools_pv"
db_name <- "refdb_teleo_2023-07-31"
fasta_name <- "input.fasta"
taxonomy_name <- "taxdump_2023-07-24.tar.gz"
db_namepath <- paste0(db_name, "/", db_name) 
#db_namepath <- "test/test2" #enkel voor de test

e <- 5
l <- 20
L <- 100
p1 <- select_primer_tags("TELEO")["oligo1"]
p2 <- select_primer_tags("TELEO")['oligo2']


## Maak structuur aan en kopieer input.fasta
system2("docker", args = c("exec", docker_container, "mkdir", db_name))
system2("docker", args = c("cp", file.path("database", db_name, fasta_name), 
                           paste0(docker_container,":/",db_name)))

#KOPIEREN TAXONOMY indien nog niet gebeurd
system2("docker", args = c("cp", file.path("taxonomy", taxonomy_name), 
                           paste0(docker_container,":/", "taxdump.tar.gz")))
##LET OP/ INDIEN NOG NIET GEBEURD MOET DE TAXONOMY UITGEPAKT WORDEN in docker met tar -xzvf taxdump.tar.gz

## converteer input+taxonomie naar voor obitools bruikbaar formaat
system2("docker", 
        args = c("exec", docker_container, "obiconvert", "--fasta", 
                        paste0("--ecopcrdb-output=",db_namepath),
                        "-t", "taxdump.tar.gz", 
                        paste0(db_name, "/", fasta_name)))

## voer de in silico PCR uit
#als ecoPCR niet gekend is ./home/ecopcr/src/ecoPCR
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

#enkel voor de case met test
# system2("docker", args = c("cp", paste0(docker_container, ':', "test"),
#                            file.path("database")))

