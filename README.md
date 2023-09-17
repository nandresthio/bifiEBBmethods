# Bi-fidelity Surrogate Modelling methods


The software in this repository and data linked to it provide methods for Bi-Fidelity Expensive Black Box (Bf-EBB) problems, both for surrogate model building and for optimisation. It implements both Kriging and Co-Kriging, two classical surrogate model methods from the literature. It relies on the benchmark repository [bifiEBBbenchmarks](https://github.com/nandresthio/bifiEBBbenchmarks) to define bi-fidelity functions, and assesses the performance of different techniques in a variety of Bf-EBB problems. 


## Cite

To cite this software, please cite the paper (ONCE IT HAS BEEN PUBLISHED) using its DOI. To cite the software itself, use the following DOI.

[![DOI](https://zenodo.org/badge/661430752.svg)](https://zenodo.org/badge/latestdoi/661430752)

Below is the BibTex for citing this version of the code. (UPDATE ZENODO)

```
@misc{nandres2022,
  author = {Andr\'{e}s-Thi\'{o}, Nicolau},
  title = {Bifidelity Surrogate Modelling methods},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  note = {available for download at https://github.com/nandresthio/bifiEBBmethods},
  doi = {10.5281/zenodo.8353700}
} 
```

## Description

This code implements methods for Bf-EBB problems. Currently the code supports only the creation of Kriging and Co-Kriging models given a sample of a bi-fidelity source. Further work will include a variety of acquisition functions, different auxiliary algorithms, and new surrogate modelling techniques. Once the performance of the chosen surrogate modelling techniques has been assessed, the R code in this repository performs an analysis of the data using Instance Space Analysis[^1] as well as an instance filtering technique to remove bias[^2]. In particular, the R file `characterisingHarmfulDataSources.R` gives the processed followed to derive the results presented in "Characterising harmful data sources when constructing
multi-fidelity surrogate models" by Nicolau Andres-Thio, Mario Andres Munoz and Kate Smith-Miles. The process consists of comparing the performance of Kriging and Co-Kriging models for a variety of test problems using the Wilcoxon test, and constructing an instance space to visualise the surrogate model performance and derive new literature guidelines. The resulting constructed instance space can be seen in the folder `data/characterisingHarmfulDataSourcesISA/` available in figshare.


### Compilation

Before compiling and running the code, it is recommended to download the populated data folder available from [figshare](https://figshare.com/articles/dataset/bifiEBBmethods_data_folder/24153342). The folder structure of the data folder in this repository is made so that all the code should run without issues, but in order to replicate results and rely on already optimised sampling plans it is recommended to replace the given data folder with the one available in figshare. 

The code relies on the repository bifiEBBbenchmarks. Before compiling the C++ code, this repository must be added to the folder by following the following steps:
  - From the root directory, call `cd cpp_code` to move to this folder
  - Call `git clone https://github.com/nandresthio/bifiEBBbenchmarks.git` to get the code

One of the instances available relies on the SOLAR simulation engine, the code of which is not provided here. If these type of instances are required, add the SOLAR code following these steps:
  
  - From the root directory, call `cd cpp_code` to move to this folder
  - Call `git clone https://github.com/bbopt/solar.git` to get the SOLAR code
  - Call `cd solar/src/`
  - Call `make` 
  - Call `cd ..`
  - Call `bin\solar -check`

If all goes well, SOLAR should be installed and the executable ready to be called from the software provided here.

To compile the C++ code given here, first make sure the library [Eigen has been installed](https://eigen.tuxfamily.org/index.php?title=Main_Page). Then go to the relevant makefile (i.e. `cpp_code\Makefile.linux` or `cpp_code\Makefile.windows`) and edit the `CFLAGS` to point to the right folder (which contains the installed Eigen library). To compile the C++ code, simply move to the `cpp_code` directory and call `make linux` if using Linux, or `make windows` if using windows. Note all code (both the R files using Rstudio or Rscript, or calling the c++ executable) should be called from the root folder due to the folder structure.


[^1]: Smith-Miles, Kate, and Mario Andrés Muñoz. "Instance space analysis for algorithm testing: Methodology and software tools." ACM Computing Surveys 55.12 (2023): 1-31.
[^2]: Alipour, Hossein, Mario Andrés Muñoz, and Kate Smith-Miles. "Enhanced instance space analysis for the maximum flow problem." European Journal of Operational Research 304.2 (2023): 411-428.