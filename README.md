<p align="center">
  <img src="docs/image/logo.png" alt="MAGE Logo" width="150"/>
</p>

# MAGE (MAtlab GEnetics)

**MAGE** is a comprehensive MATLAB toolbox designed to streamline genetic data analysis workflows. It provides user‐friendly wrappers around popular genetics software, seamless integration with R and Python tools, and built-in functions for data handling, analysis, and visualization—enabling you to run end-to-end pipelines entirely from within MATLAB.

---

## Table of Contents

* [Overview](#overview)
* [Key Features](#key-features)
* [Installation](#installation)
* [Dependencies](#dependencies)
* [Documentation](#documentation)
* [License](#license)

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

## Key Features

* **GWAS Pipeline (REGENIE)**

  * Wrapper functions that invoke REGENIE directly from MATLAB
  * Automated pre-processing (e.g., sample/variant filters) and result parsing

* **Fine-Mapping & Colocalization**

  * PolyFun and susie integration for credible set construction
  * Coloc wrapper (R) for coloc-style posterior probability estimates

* **Variant Annotation**

  * Built-in functions to call Ensembl VEP and parse consequences/annotations

* **Public API Access**

  * Query gnomAD allele frequencies and constraint metrics
  * Retrieve gene/transcript information from Ensembl REST

* **Enrichment Analysis**

  * GSEApy-based GSEA and Over-Representation Analysis (ORA) through Python bridge

* **Meta-Analysis**

  * METAL wrapper for fixed/random effects meta-analysis across cohorts

* **Conditional Analysis**

  * COJO wrapper for estimating conditional and joint effects

* **Data Import**

  * Functions to read BGEN or PLINK files directly into MATLAB tables or structs
  * Export formatted output for R survival analyses

* **Performance**

  * Optional GPU acceleration for matrix operations
  * Parallel (CPU) support using MATLAB’s Parallel Computing Toolbox

* **Plotting Utilities**

  * Regional association plots (zoomed‐in view around a locus)
  * Manhattan plots with interactive zoom and annotation
  * Effect size (forest) plots for meta-analysis results

* **UK Biobank Utilities**

  * Helper functions to prepare phenotypes (continuous, case/control, time-to-event) from UKBB sample data

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

   * [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
   * [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) (optional)
   * [Database Toolbox](https://www.mathworks.com/products/database.html) (only if direct SQL/NoSQL queries are needed)

4. **Install external dependencies:**

   * **REGENIE** (for GWAS): download and install from [REGENIE releases](https://github.com/rgcgithub/regenie/releases)
   * **METAL** (for meta-analysis): download from [METAL website](http://csg.sph.umich.edu/abecasis/Metal/)
   * **GCTA/COJO** (for conditional analysis): download [GCTA binaries](http://cnsgenomics.com/software/gcta/)
   * **Python (3.7+)** with [GSEApy](https://gseapy.readthedocs.io):

     ```bash
     pip install gseapy
     ```
   * **R (4.0+)** with `coloc` package:

     ```r
     install.packages("coloc")
     ```
   * **Ensembl VEP**: follow instructions on [VEP documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)

5. **Set up environment variables (optional but recommended):**

   * Add paths to REGENIE, METAL, GCTA, and VEP executables in your system’s `PATH` or specify them in `mage_config.m`.

---

## Dependencies

* **MATLAB (R2021a or later)**
* **Required MATLAB Toolboxes**

  * Statistics and Machine Learning Toolbox
  * (Optional) Parallel Computing Toolbox
* **External Software & Packages**

  * REGENIE (GWAS)
  * METAL (Meta-analysis)
  * GCTA/COJO (Conditional analysis)
  * Ensembl VEP (Variant annotation)
  * Python (≥3.7) + GSEApy (Enrichment analysis)
  * R (≥4.0) + coloc package (Colocalization)
* **System Requirements**

  * Linux/macOS/Windows with MATLAB support
  * ≥8 GB RAM (≥16 GB recommended for large BGEN files)
  * (Optional) GPU with CUDA support for accelerated operations

---

## Documentation

Full documentation, detailed tutorials, and additional examples are available in the [**docs/**](docs/) directory and on the project’s GitHub Pages site:
[https://yourusername.github.io/MAGE/](https://yourusername.github.io/MAGE/)

---

## License

MAGE is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

<p align="center">
  &copy; 2025 Your Name | [GitHub](https://github.com/yourusername/MAGE) | [Contact](mailto:youremail@domain.com)
</p>
