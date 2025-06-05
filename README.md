<p align="center">
  <img src="docs/image/logo.png" alt="MAGE Logo" width="150"/>
</p>

# MAGE (MAtlab GEnetics)

**MAGE** is a comprehensive MATLAB toolbox designed to streamline genetic data analysis workflows. It provides user‐friendly wrappers around popular genetics software, seamless integration with R and Python tools, and built-in functions for data handling, analysis, and visualization—enabling you to run end-to-end pipelines entirely from within MATLAB.

---

## Overview

MAGE (MAtlab GEnetics) is built to make large‐scale genetic analyses accessible from a single, MATLAB‐centric interface. Rather than switching between multiple command‐line tools or forcing data exports, MAGE “wraps” standard software (e.g., REGENIE, METAL, COJO, PolyFun, susie, coloc, VEP) and provides MATLAB functions that:

* Automate the full GWAS pipeline (including REGENIE calls)
* Perform fine-mapping (PolyFun, susie) and colocalization (coloc in R)
* Annotate variants through VEP
* Query public genetics APIs (gnomAD, Ensembl REST) directly from MATLAB
* Run enrichment analyses (GSEA or Over-Representation Analysis via GSEApy)
* Conduct meta-analysis using METAL wrappers
* Carry out conditional analyses (COJO)
* Import BGEN/PLINK files directly into MATLAB for downstream analysis (GWAS, survival analysis in R)
* Leverage GPU and multi-core CPU parallelization when available
* Generate publication-ready plots (regional association, Manhattan, effect size plots)
* Simplify UK Biobank phenotype preparation for continuous, binary, or time-to-event traits

By packaging all these capabilities into one toolbox, MAGE enables you to write concise MATLAB scripts or functions that execute end-to-end workflows—without manually juggling multiple languages or file formats.

---


## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/yourusername/MAGE.git
   cd MAGE
   ```

2. **Add MAGE to your MATLAB path:**

   ```matlab
   addpath(genpath('/path/to/MAGE'))
   savepath
   ```

3. **Install required MATLAB toolboxes:**

   * Follow the instructions to get [MATLAB Package Manager](https://mathworks.com/help/install/ug/get-mpm-os-command-line.html)
   * Set the **destinationFolder** inside mpm_input_r2024b.txt (uncomment it)
   * Run `./mpm install --inputfile=/home/matlab/Documents/MATLAB/MAGE/mpm_input_r2024b.txt"`
	* On UKB-RAP `ttyd`, you can run `MAGE_RAP_INITIALIZER.m` function directly. 

---

## UK Biobank DNAnexus Research Analysis Platform

If you intend to run analyses on UK Biobank data, follow the instructions on [DNAnexus](https://community.ukbiobank.ac.uk/hc/en-gb/articles/24999460813597-Working-with-MATLAB) for setting up MATLAB on the UK Biobank Research Analysis Platform (UKB-RAP). For more details, follow the steps below within the `ttyd` instance:

- **Prepare the auth‐token JSON file for the `dx` toolkit**  
  ```bash
  export DX_TOKEN_FILE="$HOME/.dnanexus_config/environment.json"
  chmod 644 "$DX_TOKEN_FILE"
````

* **Export the project ID for Docker**

  ```bash
  export DX_PROJECT_ID="$DX_PROJECT_CONTEXT_ID"
  ```

* **Prepare the Dockerfile**

  ```bash
  cat << EOF > Dockerfile
  FROM mathworks/matlab-deep-learning:r2024b

  USER root

  # Install Python and required libraries
  RUN apt-get update && \
      apt-get install -y \
        curl \
        python3 \
        python3-pip \
        python3-dev && \
      pip3 install --upgrade pip && \
      pip3 install dxpy pandas lxml

  # Set environment variable so MATLAB finds libpython
  ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:\$LD_LIBRARY_PATH

  # Prepare DNAnexus config path
  RUN mkdir -p /home/matlab/.dnanexus_config && \
      chown -R matlab:matlab /home/matlab/.dnanexus_config

  USER matlab
  EOF
  ```

* **Build the Docker image**

  ```bash
  docker build -t matlab-dnanexus:r2024b .
  ```

* **Run the Docker container** (adjust `--gpus all` as needed for your `ttyd` instance)

  ```bash
  docker run --gpus all -it \
    -v "/mnt/project":/dnax \
    -v "$DX_TOKEN_FILE":/home/matlab/.dnanexus_config/environment.json \
    -e DX_CONFIG_FILE=/home/matlab/.dnanexus_config/environment.json \
    -e DX_WORKSPACE_ID="$DX_PROJECT_CONTEXT_ID" \
    -p 8081:6080 \
    --shm-size=64gb \
    matlab-dnanexus:r2024b -vnc
  ```

---
