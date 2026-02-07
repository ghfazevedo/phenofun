# Phenotypic character fixation probability under neutrality (PhenoFUN)

This is a python tool to simulate a hypothetical nuclear locus that controls a phenotypic character under the coalescent process and to calculate the probability of this character to be fixed in alternate states on different populations. The method was first proposed by [Masta & Maddison (2002)](https://doi.org/10.1073/pnas.072493099) and this program was used in [Azevedo et al. (in prep)](). The model assumes that the phenotypic states of the character is controlled by a single mutation in one locus and that the mutation rate is the slowest possible (parsimony).

If you use this program, please cite [Azevedo et al. (in prep)]() and refer to [this github page](https://github.com/ghfazevedo/phenofun).

## Installation

This program uses [DendroPy](https://jeetsukumaran.github.io/DendroPy) library that is installed automatically as a dependency.
Please cite [DendroPy](https://jeetsukumaran.github.io/DendroPy).

To install PhenoFUN clone this github page and use pip.

```
git clone https://github.com/ghfazevedo/phenofun
cd phenofun
pip install .
```

## Example usage
```
phenofun -t data/tree.nwck -n 1000 -s 10,10,10,10 -os 3
```
## Outputs
The program print the probability of the target *s* to the terminal and creates the files:
1. [S_statsProbs](phenofun_out/S_statsProbs.txt) with the probability of the observed s statistics provided.
2. [simulated_gene_trees.nwck](phenofun_out/simulated_gene_trees.nwck) with all simulated trees.
3. [target_gene_trees.nwck](phenofun_out/target_gene_trees.nwck) with only gene trees that show the target *s*.
4. [simulated_s.csv](phenofun_out/simulated_s.csv) with all values of *s* for all simulations.
5. [histogram.pdf](phenofun_out/histogram.pdf) and [histogram.png](phenofun_out/histogram.png) which are the histogram plots with PDF estimated curve, with the inferior 5% inferior percentile marked in red, and with a vertical red line showing the target *s* value. 
![histogram.png](phenofun_out/histogram.png)
6. [barplot.pdf](phenofun_out/histogram.pdf) and [barplot.png](phenofun_out/histogram.png) which are the plot with probability of drift and selection (non-drift) forces fixing the alternate states. 
![barplot.png](phenofun_out/barplot.png)

## Command Options

```
usage: phenofun [-h] -t TREE [-o OUT_DIR] [-n N_SIMULATIONS] -s N_SAMPLED_INDIVIDUALS -os OBSERVED_S_STATISTICS

Calculates the probability of fixation of differences in a hypothetical nuclear locus that controls phenotype under neutral divergence.

options:
  -h, --help            show this help message and exit
  -t, --tree TREE       Path to the population/species tree with population size as branch annotations (Nexus format)
  -o, --out_dir OUT_DIR
                        Output directory
  -n, --n_simulations N_SIMULATIONS
                        Number of gene trees to simulate.
  -s, --n_sampled_individuals N_SAMPLED_INDIVIDUALS
                        The number of individuals per population/species separated by comma. It should be in the same order as the populations appear in the species tree file.
  -os, --observed_s_statistics OBSERVED_S_STATISTICS
                        The target s statistics as observed in the real world data to calculate the probability of generating it through the simulations.
```

