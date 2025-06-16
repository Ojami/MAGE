###############################################################################
# Dockerfile: mage:r2024b  – lean, sign-in-friendly, no layer bloat
###############################################################################

# --------------------------------------------------------------------------- #
# 1) Base image                                                               #
# --------------------------------------------------------------------------- #
FROM mathworks/matlab-deep-learning:r2024b

# --------------------------------------------------------------------------- #
# 2) Build-time arguments                                                     #
# --------------------------------------------------------------------------- #
ARG GIT_REPO_URL
ARG GIT_BRANCH=main
ARG MLM_LICENSE_FILE        # optional, only if you license *at build time*

# --------------------------------------------------------------------------- #
# 3) System packages + Python utilities                                       #
# --------------------------------------------------------------------------- #
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        git wget curl python3 python3-pip python3-dev \
    && pip3 install --upgrade pip dxpy pandas lxml \
    && rm -rf /var/lib/apt/lists/*

# DNAnexus helper directory and runtime libraries
ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
RUN mkdir -p /home/matlab/.dnanexus_config && \
    chown -R matlab:matlab /home/matlab/.dnanexus_config

# --------------------------------------------------------------------------- #
# 4) Clone your repo with correct ownership (NO recursive chown later)        #
# --------------------------------------------------------------------------- #
USER matlab
WORKDIR /home/matlab
RUN git clone --depth=1 --branch ${GIT_BRANCH} ${GIT_REPO_URL} MAGE

# Put MAGE on the MATLAB search path (keeps exec bits exactly as in the repo)
ENV MATLABPATH=/home/matlab/MAGE

# --------------------------------------------------------------------------- #
# 5) (Optional) Install additional toolboxes with mpm (no chmod, no chown)    #
# --------------------------------------------------------------------------- #
USER root
ENV MLM_LICENSE_FILE=${MLM_LICENSE_FILE}
RUN wget -q https://mathworks.com/mpm/glnxa64/mpm -O /tmp/mpm && \
    chmod +x /tmp/mpm && \
    printf 'destinationFolder=/opt/matlab/R2024b\n' >/tmp/mpm_input && \
    cat /home/matlab/MAGE/mpm_input_r2024b.txt >> /tmp/mpm_input && \
    /tmp/mpm install --inputfile /tmp/mpm_input && \
    rm /tmp/mpm /tmp/mpm_input

# --------------------------------------------------------------------------- #
# 6) Final user and working directory                                         #
# --------------------------------------------------------------------------- #
USER matlab
WORKDIR /home/matlab
