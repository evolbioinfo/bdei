# Performance on simulated data

## Simulated data

We assessed the performance of our estimator on two data sets: (1) __medium__, a data set of 
100 medium-sized trees (200 − 500 tips) from [Voznica _et al._ 2021](https://doi.org/10.1101/2021.03.11.435006), 
and (2) __large__, a data set of 100 large trees (5 000 − 10 000 tips) and 100 large forests (5 000 − 10 000
tips).

### Medium data set
The data were downloaded from [github.com/evolbioinfo/phylodeep](https://github.com/evolbioinfo/phylodeep) and include:
 * trees: [trees.nwk](med/trees.nwk), split by tree in [tree.[0-99].nwk](large/trees)
 * real parameter values: [TARGET.csv](med/TARGET.csv), split by tree in [tree.[0-99].log](large/trees)
 * parameter values estimated by BEAST2: [BEAST2.csv](med/BEAST2.csv)

To produce these data, Voznica _et al._ generated 10 000 trees with 200 − 500 tips under BDEI model, 
with the parameter values sampled uniformly at random within the following boundaries: 
 * incubation period 1/µ ∈ [0.2, 50]
 * basic reproductive number R<sub>0</sub> = λ/ψ ∈ [1, 5]
 * infectious period ψ ∈ [1, 10]. 

They randomly selected 100 out of those 10 000 trees to evaluate them with the gold standard method, BEAST2. 
As the BDEI model requires one of the parameters to be fixed in order to become asymptomatically identifiable, 
Voznica _et al._ fixed ρ to the real value.

### Large data set

To assess the performance of our method on larger trees, we generated a data set of 100 trees of 5 000 − 10 000 tips 
(for the same parameter values as the trees in the medium data set):
 * trees: [tree.[0-99].nwk](huge/trees)

To show the applicability to forests, we have additionally generated a data set of 100 forests 
corresponding to the large trees. 
For each simulated tree in the large tree set, 
we fixed the BDEI parameters for the forest simulation to those used to simulate the tree, 
and set the simulation time to T<sub>forest</sub> = 3/4 T<sub>tree</sub>. 
We then kept adding trees simulated over T<sub>forest</sub> 
till their total number of sampled tips n<sub>total</sub> reached the minimal number of tips (5 000). 
During this process we also kept track of the number of hidden trees u 
(i.e. those that did not get any of their tips sampled over the time T<sub>forest</sub>). 
If the total number of tips exceeded the maximal number of tips (10 000) 
– i.e. before simulating the last tree there were ntotal < 5 000 sampled tips and the last
simulated tree had n<sub>tree</sub> > 10 000 − n<sub>total</sub> sampled tips over the time T<sub>forest</sub> 
– we restarted the forest simulation procedure with the T<sub>forest</sub> time reduced by 10%:
 * forests: [forest.[0-99].nwk](huge/forests)

## Data preparation pipeline 

The [Snakemake_data](Snakemake_data) contains 
a Snakemake [[Köster *et al.*, 2012](https://doi.org/10.1093/bioinformatics/bts480)] pipeline 
that splits the trees of the medium data set into separate files (one per tree), and generates the large data set.

It can be rerun as:
```bash
snakemake --snakefile Snakefile_data --keep-going 
```


