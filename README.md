# PyBDEI

Tools for fast and accurate maximum likelihood estimation
of Birth-Death Exposed-Infectious (BDEI) epidemiological
model parameters from phylogenetic trees.

The birth-death exposed-infectious (BDEI) model [[Stadler _et al._ 2014]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4205153/) 
describes the transmission of pathogens 
that feature an incubation period (when the host is already infected but not yet infectious), 
for example Ebola or SARS-CoV-2. In a phylodynamics framework, it allows to infer such epidemiological
parameters as the basic reproduction number R<sub>0</sub>, incubation period and infectious time 
from a phylogenetic tree (a genealogy of pathogen sequences). 

This implementation of the BDEI model solves the computational bottlenecks (due to high complexity of differential equations used in phylodynamics models,
previous implementations [[Stadler and Bonhoeffer 2013](https://royalsocietypublishing.org/doi/10.1098/rstb.2012.0198) and [Barido-Sottani _et al._ 2018](https://doi.org/10.1101/440982) ] sometimes suffered from numerical instability and were only applicable to medium datasets of <500 samples). 
Our fast and accurate estimator is applicable to very large datasets (10, 000 samples) allowing phylodynamics to
catch up with pathogen sequencing efforts.

#### Article

A Zhukova, F Hecht, Y Maday, and O Gascuel. *Fast and accurate maximum likelihood estimation
of Birth-Death Exposed-Infectious (BDEI) epidemiological
model from phylogenetic trees* _(in prep.)_

# Input data
As an input, one needs to provide a **rooted** phylogenetical tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format,
and the value of one of the model parameters (for identifiability):
* µ – becoming infectious rate corresponding to a state transition from E (exposed) to I (infectious) 
_(can be fixed via the --mu argument)_,
* λ – transmission rate, from a transmitter in the state I to a newly infected recipient, whose state is E 
_(can be fixed via the --la argument)_,
* ψ – removal rate, corresponding to individuals in the state I exiting the study 
(e.g. due to healing, death or starting a treatment) _(can be fixed via the --psi argument)_,
* ρ – sampling probability (upon removal) _(can be fixed via the --p argument)_.


# Installation

There are 4 alternative ways to run __PyBDEI__ on your computer: 
with [docker](https://www.docker.com/community-edition), 
[singularity](https://www.sylabs.io/singularity),
in python3 (only on linux systems), or via command line (only on linux systems, requires installation with python3).

## Run with docker

### Basic usage
Once [docker](https://www.docker.com/community-edition) is installed, run the following command 
(here we assume that the sampling probability value is known and fixed to 0.3):

```bash
docker run -v <path_to_the_folder_containing_the_tree>:/data:rw -t evolbioinfo/pybdei --nwk /data/<tree_file.nwk> --p 0.3 --CI_repetitions 100 --log <file_to_store_the_estimated_parameters.tab>
```

This will produce a file <file_to_store_the_estimated_parameters.tab> in the <path_to_the_folder_containing_the_tree> folder,
 containing a tab-separated table with the estimated parameter values and their CIs (can be viewed with a text editor, Excel or Libre Office Calc).

#### Help

To see advanced options, run
```bash
docker run -t evolbioinfo/pybdei -h
```

## Run with singularity

### Basic usage
Once [singularity](https://www.sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps) is installed, 
run the following command  
(here we assume that the sampling probability value is known and fixed to 0.3):

```bash
singularity run docker://evolbioinfo/pybdei --nwk <path/to/tree_file.nwk> --p 0.3 --CI_repetitions 100 --log <path/to/file_to_store_the_estimated_parameters.tab>
```

This will produce a file <path/to/file_to_store_the_estimated_parameters.tab>,
 containing a tab-separated table with the estimated parameter values and their CIs (can be viewed with a text editor, Excel or Libre Office Calc).


#### Help

To see advanced options, run
```bash
singularity run docker://evolbioinfo/pybdei -h
```

## Run in python3 or command-line (for linux systems, recommended Ubuntu 21 or older versions)

### 1. Install the C++ dependencies
You would need to install g++10 and libnlopt-cxx-dev libraries:

```bash
sudo apt update --fix-missing 
sudo apt install -y g++-10 libnlopt-cxx-dev
```

### 2. Install python 3

You could either install python 3 system-wide:
```bash
sudo apt install -y python3 python3-pip python3-setuptools python3-distutils
```

or alternatively, you could install python 3 via [conda](https://conda.io/docs/) (make sure that conda is installed first). 
Here we will create a conda environment called _pybdeienv_:
```bash
conda create --name pybdeienv python=3
conda activate pybdeienv
pip3 install setuptools
```

### 3. Install numpy and PyBDEI
```bash
pip3 install numpy 
pip3 install pybdei
```


### Basic usage in a command line
If you installed __PyBDEI__ via conda, do not forget to first activate the dedicated environment (here named _pybdeienv_), e.g.

```bash
conda activate pybdeienv
```

To run __PyBDEI__
(here we assume that the sampling probability value is known and fixed to 0.3):

```bash
bdei_infer --nwk <path/to/tree_file.nwk> --p 0.3 --CI_repetitions 100 --log <path/to/file_to_store_the_estimated_parameters.tab>
```

This will produce a file <path/to/file_to_store_the_estimated_parameters.tab>,
 containing a tab-separated table with the estimated parameter values and their CIs (can be viewed with a text editor, Excel or Libre Office Calc).

#### Help

To see advanced options, run:
```bash
bdei_infer -h
```

### Basic usage in python3

```python
from pybdei import infer
# Path to the tree in newick format
tree = "tree.nwk"
result, time = infer(nwk=tree, p=0.3, CI_repetitions=100)
print('Inferred transition rate is', result.mu, result.mu_CI)
print('Inferred transmission rate is', result.la, result.la_CI)
print('Inferred removal rate is', result.psi, result.psi_CI)
print('Inferred reproductive number is', result.R_naught)
print('Inferred incubation period is', result.incubation_period)
print('Inferred infectious time is', result.infectious_time)
print('Converged in', time.CPU_time, 's and', time.iterations, 'iterations')
```