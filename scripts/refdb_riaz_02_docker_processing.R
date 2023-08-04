#Build docker image
#Zorg dat docker desktop voor windows aan het draaien is (of het docker process in linux)

#Bouw de image ("docker build -t obitools3pv ./docker_image_obitools3/")
system2("powershell", args = pscommand_build)

#Run de container van de imag
  #delete eerst docker rm obitools3testcontainer
  #run: "docker run --rm -i -d --name obitools3testcontainer obitools3pv")
system2("powershell", args = pscommand_delete)
system2("powershell", args = pscommand_run)

#als het bestand nog niet bestaat: op de juiste locatie zetten en zorgen dat het naar docker gekopieerd wordt
#-------------------
#Maak het shellscript aan (generate_refdb.sh) die de refdb zal aanmaken
# make_shellscript_refdb(script = "generate_refdb.sh",
#                        db_location = file.path("database",db_name), 
#                        db_name = db_name,
#                        input_file = "input.fasta",
#                        taxonomy_location = file.path("taxonomy"),
#                        taxonomy_file = "taxdump.tar.gz",
#                        max_errors = 5, 
#                        min_length = 20, 
#                        max_length = 100,
#                        primer = "RIAZ",
#                        ecotag_min_similarity = 0.96)

##copy to container
system2('docker',  c('cp', file.path("docker_obitools3_Teleo", obi_script),
                     paste0(docker_container_name, ':', docker_script)))
system2('docker',  c('cp', taxdump_name, docker_taxonomy))
system2('docker',  c('cp', file.path(refdb_location, 'input.fasta'), docker_input_fasta))

##execute main script
system2('docker', c('exec', docker_container_name, 'bash', '-c', docker_script))

##return results
system2('docker', c("cp", paste0(docker_path, "logfile.txt"), refdb_location))
system2('docker', c("cp", paste0(docker_path, obi_script), refdb_location))
system2('docker', c("exec", docker_container_name, "chmod", "777", obi_script))
system2('docker', c("cp", paste0(docker_path, "kept_input.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "amplified_clean.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "amplified_clean_uniq.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "final_db_0.96.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "final_db_0.99.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "final_db_1.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "refdb.obidms"), refdb_location))        

##exit container
system2("powershell", paste("docker stop ", docker_container_name))       

