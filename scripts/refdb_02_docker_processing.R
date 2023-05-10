#Build docker image
#Zorg dat docker desktop voor windows aan het draaien is (of het docker process in linux)

#Bouw de image ("docker build -t obitools3pv ./docker_image_obitools3/")
system2("powershell", args = pscommand_build)

#Run de container van de imag
  #delete eerst docker rm obitools3testcontainer
  #run: "docker run --rm -i -d --name obitools3testcontainer obitools3pv")
system2("powershell", args = pscommand_delete)
system2("powershell", args = pscommand_run)

#Maak het shellscript aan (generate_refdb.sh) die de refdb zal aanmaken
make_shellscript_refdb(script = "generate_refdb.sh",
                       db_location = file.path("database",db_name), 
                       db_name = db_name,
                       input_file = "input.fasta",
                       taxonomy_location = file.path("taxonomy"),
                       taxonomy_file = "taxdump.tar.gz",
                       max_errors = 4, 
                       min_length = 50, 
                       max_length = 160,
                       primer = "RIAZ",
                       ecotag_min_similarity = 0.99)

##copy to container
system2('docker',  c('cp', obi_script, paste0(container_name, ':', docker_script)))
system2('docker',  c('cp', taxdump_name, docker_taxonomy))
system2('docker',  c('cp', input_fasta, docker_input_fasta))

##execute main script
system2('docker', c('exec', container_name, 'bash', '-c', docker_script))

##return results
system2('docker', c("cp", paste0(docker_path, "logfile.txt"), refdb_location))
system2('docker', c("cp", paste0(docker_path, obi_script), refdb_location))
system2('docker', c("exec", container_name, "chmod", "777", obi_script))
system2('docker', c("cp", paste0(docker_path, "kept_input.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "amplified_clean.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "amplified_clean_uniq.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "final_db_0.99.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "refdb.obidms"), refdb_location))        

##exit container
system2("powershell", paste("docker stop ", container_name))       


df_inputs %>% filter(taxid == 96503) %>% view()
