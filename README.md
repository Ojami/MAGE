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
	* On UKB-RAP `ttyd`, you can run MAGE_RAP_INITIALIZER.m function directly. 

---