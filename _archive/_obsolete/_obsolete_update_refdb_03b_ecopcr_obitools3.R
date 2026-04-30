
source("scripts/update_refdb_00_init.R")

#Obitools3:
# sudo apt-get update -y
# sudo apt install -y python3-venv
# sudo apt-get install python3-pip
# sudo apt-get install cmake
# python3 -m venv obi3-env
# source obi3-env/bin/activate
# pip3 install --upgrade pip setuptools wheel Cython
# pip3 install OBITools3

make_shellscript_refdb(script = "console",
                       db_location = "database/refdb_2022-11-11", 
                       db_name = "refdb_2022-11-11",
                       input_file = "input_obi3.fasta",
                       taxonomy_location = "taxonomy/taxdump_2022-11-09",
                       taxonomy_file = "new_taxdump_2022-11-09.tar.gz")


