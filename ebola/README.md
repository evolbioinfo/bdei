# Ebola SLE 2014 epidemic analysis

## Data

We took the 1,610 sequence alignment and metadata (sampling times and countries) 
that were used in the study by [Dudas _et al._ 2017](https://www.nature.com/articles/nature22040), 
downloaded from [github.com/ebov/space-time](https://github.com/ebov/space-time/tree/master/Data) 
and saved here as:
 * [data/aln.ids.fa](data/aln.ids.fa) -- multiple sequence alignment;
 * [data/metadata.csv](data/metadata.csv) -- metadata.
 
The original data are available under [Creative Commons Attribution Share Alike 4.0 International licence](LICENCE.txt), the same licance applies to the results that we obtained from it (see below).

### SLE forests

We reconstructed the maximum likelihood phylogeny for this alignment with RAxML-NG (v1.0.2, GTR+G4+FO+IO) [[Kozlov _et al._ 2019](https://pubmed.ncbi.nlm.nih.gov/31070718/)], 
and rooted it based on sampling dates using LSD2 (v1.10) [[To _et al._ 2016](https://academic.oup.com/sysbio/article/65/1/82/2461506)]:
 * [data/tree.nwk](data/tree.nwk) -- phylogeny
 * [data/timetree.nexus](data/timetree.nexus) -- rooted (and time-scaled) tree

As Ebola's mutation rate is slower than its transmission rate, the initial phylogeny contained 246 polytomies. BDEI
model, on the other hand, assumes a binary tree. We therefore resolved these polytomies randomly (10 times, to check for
robustness of the estimates) using a coalescent approach:
 * [data/rooted.tree.0-9.nwk](data)

We dated each of the 10 trees with LSD2 (v1.10 [[github.com/tothuhien/lsd2](github.com/tothuhien/lsd2/tree/v1.10)], 
under strict molecular clock with outlier removal) using tip sampling dates, and reconstructed the ancestral characters for country with PastML
(v1.9.34, MPPA+F81) [[Ishikawa, Zhukova _et al._ 2019](https://academic.oup.com/mbe/article/36/9/2069/5498561)]:
 * [data/timetree.0-9.nexus](data) -- time-scaled trees
 * [data/timetree.0-9.itol](data) -- iTOL [[Letunic & Bork 2019](https://academic.oup.com/nar/article/47/W1/W256/5424068?login=true)] visualisations of the ancestral country reconstructions.

Lastly, we extracted 10 SLE forests from these trees to represent the Ebola epidemic in SLE between July 30 2014 (when the SLE
government began to deploy troops to enforce quarantines according to [news24.com](https://web.archive.org/web/20190505224120/https://www.news24.com/Africa/News/Sierra-Leone-Liberia-deploy-troops-for-Ebola-20140804)) 
 and September 12 2015 (the last SLE sample in these data) by 
1. cutting each tree on July 30 2014 to remove the more ancient part (with a different health policy); 
2. among the July-31-on trees, picking those whose root’s predicted character state for country was SLE; 
3. removing the non-SLE subtrees from the selected July-31-on SLE trees to focus on the epidemic within the country, 
without further reintroductions.

For each of the forests we converted its tree branch lengths to days:
 * [data/SLE/SLE.0-9.days.nwk](data/SLE)

### Pipeline


The [Snakemake_data](Snakemake_data) file contains 
a Snakemake [[Köster *et al.*, 2012](https://doi.org/10.1093/bioinformatics/bts480)] pipeline 
that reconstructs the 10 SLE forests as described above.

It can be rerun as:
```bash
snakemake --snakefile Snakefile_data --keep-going  --config folder=data --use-singularity --singularity-prefix ~/.singularity --singularity-args "--home ~"
```

# BDEI parameter estimation

We estimated the BDEI parameters on these 10 forests, fixing the total number of trees to N = 533 (the declared number of cases in
SLE on July 31, 2014 from [cdc.gov](https://www.cdc.gov/vhf/ebola/history/2014-2016-outbreak/case-counts.html)), 
and hence the number of unobserved trees to N − k, where k was the number of trees in the corresponding
forest (k varied between 55 and 70). 
To check the robustness of the predictions with respect to N, 
we additionally estimated the parameter values assuming 50% more cases, i.e. N = 800.

As BDEI model requires one of the parameters to be fixed in order to become asymptomatically identifiable, 
we performed the estimations fixing the infectious period (to 2.6 and to 5 days, i.e. the estimates from the previous studies). 

The [Snakemake_estimate](Snakemake_estimate) file contains 
a Snakemake pipeline for parameter estimation.

It can be rerun as:
```bash
snakemake --snakefile Snakefile_estimate --keep-going  --config folder=.
```

The results are listed in [data/SLE/estimates.tab](data/SLE/estimates.tab).
