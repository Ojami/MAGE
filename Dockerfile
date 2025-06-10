###############################################################################
# Dockerfile: mage-full:r2024b (dynamic mpm_input edit)
###############################################################################

##########################################
# 1) Base Image
##########################################
FROM mathworks/matlab-deep-learning:r2024b

##########################################
# 2) Build Arguments
##########################################
ARG GIT_REPO_URL
ARG GIT_BRANCH=main
ARG MLM_LICENSE_FILE

##########################################
# 3) Install system dependencies (as root)
##########################################
USER root
RUN apt-get update && \
    apt-get install -y \
      git \
      wget \
      curl \
      python3 \
      python3-pip \
      python3-dev && \
    pip3 install --upgrade pip && \
    pip3 install dxpy pandas lxml && \
    rm -rf /var/lib/apt/lists/*

ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
RUN mkdir -p /home/matlab/.dnanexus_config && \
    chown -R matlab:matlab /home/matlab/.dnanexus_config

##########################################
# 4) Switch to matlab user, clone MAGE, download mpm
##########################################
USER matlab
WORKDIR /home/matlab

RUN git clone --branch ${GIT_BRANCH} ${GIT_REPO_URL} /home/matlab/MAGE

RUN wget https://mathworks.com/mpm/glnxa64/mpm -O /home/matlab/mpm && \
    chmod +x /home/matlab/mpm

##########################################
# 5) Switch back to root, set license, run mpm install
##########################################
USER root
ENV MLM_LICENSE_FILE=${MLM_LICENSE_FILE}

# (5a) Build a temporary mpm input file that prepends destinationFolder
RUN printf "destinationFolder=/opt/matlab/R2024b\n" > /home/matlab/tmp_mpm_input.txt && \
    cat /home/matlab/MAGE/mpm_input_r2024b.txt >> /home/matlab/tmp_mpm_input.txt

# (5b) Add MAGE to MATLAB path, then install toolboxes using that temp file
RUN matlab -batch "addpath(genpath('/home/matlab/MAGE')); savepath; exit" && \
    /home/matlab/mpm install --inputfile /home/matlab/tmp_mpm_input.txt

##########################################
# 6) Final USER and WORKDIR
##########################################
USER matlab
WORKDIR /home/matlab
CMD ["bash"]
