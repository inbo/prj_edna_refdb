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
    
    
# Update PATH variable
ENV PATH="/app/obitools3/bin:${PATH}"

# Configure OBITools3 auto-completion (optional)
RUN echo "source /app/obitools3/obi_completion_script.bash" >> ~/.bashrc

# Add obi command to PATH
RUN echo 'export PATH="$PATH:/app/obitools3/bin"' >> ~/.bashrc

# Copy your Bash script into the container
COPY welcome.sh /app/welcome.sh

#Copy to run non-interactively en execute when building the image
#COPY generate_refdb.sh /app/generate_refdb.sh
#COPY input.fasta /app/input.fasta
#COPY new_taxdump_2022-05-02.tar.gz /app/taxdump.tar.gz

# Set the command to run your Bash script
#CMD ["bash", "/app/welcome.sh"] #docker sluit na afronden script
#CMD ["/bin/bash", "-c", "/app/welcome.sh && exec /bin/bash"] #docker image blijft open
CMD ["/bin/bash", "-c", "source /app/obi3-env/bin/activate && obi test && /app/welcome.sh && exec /bin/bash"]

# Set the command to keep the container running indefinitely
CMD ["tail", "-f", "/dev/null"]


