# Phenotypic character fixation probability under neutrality (PhenoFUN)

This is a python tool to simulate a hypothetical nuclear locus that controls a phenotypic character under the coalescent process and to calculate the probability of this character to be fixed in alternate states on different populations. The method was first proposed by [Masta & Maddison (2002)](https://doi.org/10.1073/pnas.072493099) and this program was used in [Azevedo et al. (in prep)](). The model assumes that the phenotypic states of the character is controlled by a single mutation in one locus and that the mutation rate is the smallest possible (parsimony).

If you use this program, please cite [Azevedo et al. (in prep)]() and refer to github page.

## Installation

This program requires [DendroPy](https://jeetsukumaran.github.io/DendroPy).
```
pip install dendropy==5.0.1
```
After installing dependencies, clone the repository and install it

```
git clone https://github.com/ghfazevedo/phenofun
cd phenofun
pip install .
```

## Example usage
```
phenofun -t data/tree.nwck -n 100 -s 10,10,10,10 -os 3
```

## Command Options

```
usage: phenofun.py [-h] -t TREE [-o OUT_DIR] [-n N_SIMULATIONS] -s N_SAMPLED_INDIVIDUALS -os OBSERVED_S_STATISTICS

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

