# Use the official Python 3.10 image as the base image
FROM python:3.10-slim

# Set the working directory
WORKDIR /app

# Update package index and install required packages
RUN apt-get update -y && \
    apt-get install -y vim wget gcc zlib1g-dev make git

RUN apt-get update \
  && apt-get -y install build-essential \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget https://github.com/Kitware/CMake/releases/download/v3.24.1/cmake-3.24.1-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /opt/cmake-3.24.1 \
      && /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake-3.24.1 \
      && rm /tmp/cmake-install.sh \
      && ln -s /opt/cmake-3.24.1/bin/* /usr/local/bin



# Clone OBITools3 repository
RUN git clone https://git.metabarcoding.org/obitools/obitools3.git

# Create and activate virtual environment
RUN python3 -m venv obi3-env && \
    . obi3-env/bin/activate

# Install Cython
RUN obi3-env/bin/pip install cython

# Install OBITools3
RUN . obi3-env/bin/activate && \
    cd obitools3 && \
    python3 setup.py install && \
    cd .. \
    mkdir refdb

# Install Cutadapt
RUN obi3-env/bin/pip install cutadapt

#Install Flash
RUN wget https://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz && \
    tar -xzf FLASH-1.2.11.tar.gz && \
    cd FLASH-1.2.11 && \
    make && \
    cp flash /usr/local/bin/ && \
    cd .. && \
    rm -rf FLASH-1.2.11 && \
    rm FLASH-1.2.11.tar.gz

#gewoon om cutadapt te testen (R1 5'--> 3', dus behou alles na -g of alles voor -a)
#python3 -m venv obi3-env && . obi3-env/bin/activate
#cutadapt -g GTCAGGGGAAAA -o outg.fastq testfwd.fastq
#cutadapt -a GTCAGGGGAAAA -o outa.fastq testfwd.fastq
#hieronder lijkt nog niet te lukken
#cutadapt -g GTCAGGGGAAAA -a AATAAG -o outga.fastq testfwd.fastq
#cutadapt -g GTCAGGGGAAAA -a CTTATT -o outag.fastq testfwd.fastq

#cutadapt -g TAAGAAT -G CGCGCG -o outR1.fastq -p outR2.fastq testfwd.fastq testrev.fastq
#flash -m 4 outR1.fastq outR2.fastq
    
    
# Update PATH variable
ENV PATH="/app/obitools3/bin:${PATH}"

# Configure OBITools3 auto-completion (optional)
RUN echo "source /app/obitools3/obi_completion_script.bash" >> ~/.bashrc

# Add obi command to PATH
RUN echo 'export PATH="$PATH:/app/obitools3/bin"' >> ~/.bashrc

# Copy your Bash script into the container
COPY welcome.sh /app/welcome.sh

# Copy testfastq into container
#COPY testfwd.fastq /app/testfwd.fastq
#COPY testrev.fastq /app/testrev.fastq

#Copy to run non-interactively en execute when building the image
#COPY generate_refdb.sh /app/generate_refdb.sh
#COPY input.fasta /app/input.fasta
#COPY new_taxdump_2022-05-02.tar.gz /app/taxdump.tar.gz

# Set the command to run your Bash script
#CMD ["bash", "/app/welcome.sh"] #docker sluit na afronden script
#CMD ["/bin/bash", "-c", "/app/welcome.sh && exec /bin/bash"] #docker image blijft open
CMD ["/bin/bash", "-c", "source /app/obi3-env/bin/activate && /app/welcome.sh && exec /bin/bash"]

# Set the command to keep the container running indefinitely
CMD ["tail", "-f", "/dev/null"]


