# Analysis_Of_Ribosome_Pausing
A repository that details an end-to-end pipeline for performing a differential analysis of pausing in ribosome profiling data.

# Installation

## Plastid and python environment preparations.

**Timing: 20 minutes**

Plastid is a Python library which contains a variety of tools for genomics and sequencing and is available on all operating systems (Windows, Linux, Mac OSX). This library will be used extensively within this analysis. The Plastid Python library is available as an Anaconda package library and can be easily installed within a Conda environment using Miniconda. Other necessary Python packages can be subsequently installed within the Plastid environment. 
		 	 	 		
1. Install Miniconda and setup Bioconda on your system if you have not already. 

2. Create a conda environment with plastid installed:
`> conda create -n plastid plastid`

3. Enter the new conda environment:
`> conda activate plastid`

4. Plastid is designed to work with an older version of biopython and is incompatible with newer versions. Therefore, it is necessary to downgrade Biopython by running:
`> conda install biopython==1.76` 

5. Install Jupyterlab within the Plastid environment:
`> conda install -c conda-forge jupyterlab`

6. Install the multiprocess Python package using conda:
`> conda install multiprocess`

## R environment preparations. 

**Timing: 30 minutes**

the Determining P-site offsets section of this protocol requires the use of an R-package called riboWaltz. This section must be installed within R directly from GitHub using the devtools suite of R-packages. 

1. Install R and R studio using the instructions from this link https://www.r-project.org/ and this link https://www.rstudio.com/products/rstudio/download/ respectively.

2. Open R-studio and install the devtools suite of packages using R’s install.packages function and then load devtools into the R session using the library command:
 
`> install.packages("usethis")`
`> install.packages("devtools")` 
`> library(devtools)`

3. Install riboWaltz directly from GitHub using devtools' install_github function: 

  `> install_github("LabTranslationalArchitectomics/riboWaltz",
  dependencies = TRUE, build_opts = c("--no-resave-data", "-
  no-manual"))`
