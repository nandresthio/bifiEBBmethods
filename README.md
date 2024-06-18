# Bi-fidelity Surrogate Modelling methods


The software in this repository and data linked to it provide methods for Bi-Fidelity Expensive Black Box (Bf-EBB) problems, both for surrogate model building and for optimisation. It implements both Kriging and Co-Kriging, two classical surrogate model methods from the literature, which are used in the paper [Characterising Harmful Data Sources When Constructing Multi-fidelity Surrogate Models](https://arxiv.org/pdf/2403.08118) by Andrés-Thió N, Muñoz MA, and Smith-Miles K to characterise when additional low-fidelity sources are benefitial and when they are detrimental in surrogate modelling construction. R and MATLAB code is also provided to process the performance results of these techniques.

Furthermore, implementations are provided for the acquisition function known as Integrated Mean Square Prediction Error (IMSPE) for both Kriging and Co-Kriging, the theory for which has been newly developed and is presented in the paper [TO BE DEFINED](link.com). In the same paper defines the methods Rules-based Co-Kriging, Adaptive Co-Kriging, and Rules-based Adaptive Co-Kriging, all of which are implemented here.

Note that this repository relies on the benchmark repository [bifiEBBbenchmarks](https://github.com/nandresthio/bifiEBBbenchmarks) to define bi-fidelity functions; further details on the compilation are given below.


## Cite

This software can be of interest for a variety of reasons. If it is used to assess surrogate models given a fixed sample (by rerunning the code or using the data available on figshare), following the experimental setup of the paper [Characterising Harmful Data Sources When Constructing Multi-fidelity Surrogate Models](https://arxiv.org/pdf/2403.08118), please cite the preprint paper using the following BibTex citation:

```
@article{andres2024characterising,
  title={Characterising harmful data sources when constructing multi-fidelity surrogate models},
  author={Andr{\'e}s-Thi{\'o}, Nicolau and Mu{\~n}oz, Mario Andr{\'e}s and Smith-Miles, Kate},
  journal={arXiv preprint arXiv:2403.08118},
  year={2024}
}
```

If this software is used due to the implemented techniques proposed in the paper [TO BE DEFINED](link.com) such as the acquisition function IMSPE, please cite the preprint paper using the following BibTex citation:

ADD WHEN AVAILABLE

To cite the software itself, use the following DOI: [![DOI](https://zenodo.org/badge/661430752.svg)](https://zenodo.org/badge/latestdoi/661430752); below is the BibTex for citing this version of the code.

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

This code implements methods for Bi-Fidelity Expensive Black-Box problems. These problems have a high- and low-fidelity source of information denoted $f_h$ and $f_l$ respectively, which are used to either accurately model or to optimise the high-fidelity source. The code provides an implementation of the classical surrogate modelling techniques known as Kriging and Co-Kriging, as well as newly proposed techniques which are described below. The implementation also provides code which reads in an experiment specification (also detailed below) and assesses the performance of the specified technique on the specified problem. 


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


### Available techniques

Surrogate modelling techniques consist of methods which train a model of a given black-box $f_h$. This model may be used to guide further sampling, but it is possible for the construction of an accurate model to be the end-goal. This repository implements the classical techniques Kriging[^1] and Co-Kriging[^2]. Kriging is a single-source model, meaning it trains a predictive model of $f_h$ using data from this source exclusively. Co-Kriging is a two-source technique, meaning it incorporates data from $f_h$ as well as a low-fidelity source $f_l$ so generate predictions of $f_h$. Both of these techniques are described in detail in the thesis provided in `thesis` folder, in Appendix A.

The same thesis proposes Rules-based Co-Kriging. This surrogate modelling approach is simply an additional layer which assesses whether a Kriging model or a Co-Kriging model should be trained, given the characteristics of the available data. Details on this approach can be found in the provided thesis in Chatper 5.

Given a trained surrogate model, further sampling (either to increase the accuracy of a model or for optimisation) is decided via an acquisition function. Classical acquisition functions for both Kriging and Co-Kriging include Expected Improvement (when optimising) and Maximum Prediction Error (MPE) (when modelling). Both of these are implemented in this repository. It is also possible to sample by minimising the overall uncertainty of future models (known as IMSPE); doing so has traditionally been computationally expensive. In this repository however and for the first time, an efficient implementation has been derived as is presented in Chapter 5 of the provided thesis.

In terms of generating an initial design of experiments, this repo implements the approach suggested by Forrester et al.[^2]. This consists of first generating an initial Latin Hypercube Sampling (LHS) plan, and to then locally optimise it so that the samples are spread out. A similar approach is taken to choose a subset of the samples to generate a nested sampling plan. Note than this approach can be lengthy for larger problem dimensions and sample sizes. Due to this and as stated above it is recommended to download the populated data folder available from [figshare](https://figshare.com/articles/dataset/bifiEBBmethods_data_folder/24153342) prior to running the code.


### Testing

This repository is primarily made available for future researchers to validate the work presented in [Characterising Harmful Data Sources When Constructing Multi-fidelity Surrogate Models](https://arxiv.org/pdf/2403.08118) and [TO BE DEFINED](link.com). As such, the performance of the surrogate modelling methods can be assessed on the following Bf-EBB subproblems, both of which are defined in Chapter 2 of the provided thesis:

 - Model creation with fixed sample: This is defined by the tuple $(f_h, f_l, n_h, n_l)$, where $f_h$ is the high-fidelity source to be modelled, $f_l$ is an available low-fidelity source of unkown quality, and $n_h$ and $n_l$ denote the number of available high- and low-fidelity sources, respectively. The chosen surrogate modelling method does not get to gather further samples, and only needs to train an accurate model given the available data. The performance can either be measured in terms of error or correlation with $f_h$; both are recorded by the code.

 - Model creation with sample budget: This is defined by the tuple $(f_h, f_l, B, C_r)$, where $f_h$ and $f_l$ are defined as above, $B$ is the total number of high-fidelity evaluations available, and $C_r$ is the cost of sampling $f_l$ relative to $f_h$. This means having sampled $f_h$ and $f_l$ a total of $n_h$ and $n_l$ times, respectively, the overall budget used is given by $n_h + C_r n_l$. The chosen technique needs to first choose how much of the total budget use on an initial design of experiments, and how to gather further data samples.

In order to run a set of tests, the compile code (compiled following the instructions given below) should be run by specifying by calling the executable, along with the name of the file containing the test run information and the line in the file to be run. An example run may be `cpp_code\main experimentalRunSurrogateModelWithFixedSample 25 10` which specifies to run the experiments defined in `data\runScripts\experimentalRunSurrogateModelWithFixedSample.txt`, and to run line 25 followed by the next 10 runs.

The files stored in `data\runScripts` prove examples on how to define an experimental run. For the first sub problem, the first column should specify `surrogateModelWithFixedSample`. The second column should define the instance. this is defined by an instance name, followed by $n_h$ and $n_l$, followed by which seeds to run. An example is `(ToalPaciorek0.60,4,8,1-40)`. The third column should specify the technique being tested, which should be one of `kriging` or `cokriging`.

If testing the second subproblem, the first column should specify this with `surrogateModelWithBudget`. The second column defines the instance, which here is defined by an instance name, $B$, $C_r$ and the seeds to run. An example is `(ToalPaciorek0.60,20,0.1,1-40)`. Note that $B$ is given in relative terms; that is in the given example a total, as Paciorek is a 2-dimensional function, a total budget of $B = 20 \times 2 = 40$ will be used. The last column specifies the technique to use in the form of `MODEL_ACQUISITION_DOE`. Here `MODEL` is one of `kriging`, `cokriging` and `rulesCokriging` (the last of which corresponds to rules-based Co-Kriging as defined in Chapter 5 of the thesis), `ACQUISITION` defines the acquisition function and should be one of `variance` (for maximising model uncertainty), `globalVariance` (for IMSPE) or `globalVarianceWithChoice` (for Adaptive Co-Kriging as defined in Chapter 5of the thesis), and `DOE` specifies how much of the total budget to use on an initial DoE: either `all` (for all of it), `half` (for half of it), or `small` (for $n_h = 2d + 1$ where $d$ is the problem dimension). Note that in all three cases twice as many low-fidelity samples are generated as high-fidelity samples for the DoE.

The testing results are stored in the folder `data\clusterResults` to be processed as detailed below.


### Results processing

The R files stored in the folder `R_code` are used to process and analyse all the results as presented in the two papers discussed in this repo.

For the paper [Characterising Harmful Data Sources When Constructing Multi-fidelity Surrogate Models](https://arxiv.org/pdf/2403.08118), the results are processed using the Wilcoxon test as detailed both in the paper and [here](data/isaMetadata). They are then analysed using Instance Space Analysis (ISA)[^3] as well as an instance filtering technique to remove bias[^4]. Note that this repository provides code for both of these techniques and the file `R_code/characterisingHarmfulDataSources.R` gives the process followed to derive the results presented in the paper. The reader interested in these specific techniques however should go to the [original ISA repo](https://github.com/andremun/InstanceSpace) and cite accordingly. The resulting constructed instance space can be seen in the folder `data/characterisingHarmfulDataSourcesISA/` available in figshare.

For the paper [TO BE DEFINED](link.com) the processing of results is simply aggregated by budget and FINISH ONCE THIS IS FINALISED.


## Acknowledgements

This research was supported by the Australian Research Council under grant number IC200100009 for the ARC Training Centre in Optimisation Technologies, Integrated Methodologies and Applications (OPTIMA). Funding was also received through the Australian Government Research Training Program Scholarship from The University of Melbourne. This research was also supported by The University of Melbourne’s Research Computing Services and the Petascale Campus Initiative.




[^1]: Donald R Jones. A taxonomy of global optimization methods based on response surfaces. Journal of global optimization, 21(4):345–383, 2001.
[^2]: Alexander IJ Forrester, Andr´as S´obester, and Andy J Keane. Multi-fidelity optimization via surrogate modelling. Proceedings of the royal society a: mathematical,
physical and engineering sciences, 463(2088):3251–3269, 2007.
[^3]: Smith-Miles, Kate, and Mario Andrés Muñoz. "Instance space analysis for algorithm testing: Methodology and software tools." ACM Computing Surveys 55.12 (2023): 1-31.
[^4]: Alipour, Hossein, Mario Andrés Muñoz, and Kate Smith-Miles. "Enhanced instance space analysis for the maximum flow problem." European Journal of Operational Research 304.2 (2023): 411-428.