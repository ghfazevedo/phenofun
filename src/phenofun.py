#!/usr/bin/env python3

import argparse
import os
import dendropy
from dendropy.simulate import treesim
from dendropy.model import reconcile
from dendropy.model import coalescent
from dendropy.model.reconcile import monophyletic_partition_discordance
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from phenofun import __version__
#from scipy.stats import lognorm

# Function to parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculates the probability of fixation of differences in a hypothetical nuclear locus that controls phenotype under neutral divergence.")
    parser.add_argument("-t", "--tree", required=True, help="Path to the population/species tree with population size as branch annotations (Nexus format)")
    parser.add_argument("-o", "--out_dir", default="./phenofun_out", help="Output directory")
    parser.add_argument("-n", "--n_simulations", default="100", help="Number of gene trees to simulate.")
    parser.add_argument("-s", "--n_sampled_individuals", type=str, required=True, help="The number of individuals per population/species separated by comma. It should be in the same order as the populations appear in the species tree file.")
    parser.add_argument("-os", "--observed_s_statistics", type=int, required=True, help="The target s statistics as observed in the real world data to calculate the probability of generating it through the simulations.")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    return parser.parse_args()

# Function to confirm with the user if directory exists
def confirm_proceed(message="Directory already exists. Do you want to proceed? (y/n): "):
    while True:
        response = input(message).strip().lower()
        if response == 'y':
            return True
        elif response == 'n':
            print("Exiting program.")
            return False
        else:
            print("Please enter 'y' or 'n'.")

def main():
    args = parse_arguments()

    # Convert out_dir to absolute path
    args.tree = os.path.abspath(args.tree)
    args.out_dir = os.path.abspath(args.out_dir)

    # Convert the number of individuals to a list
    n_sampled_individuals = list(map(int, args.n_sampled_individuals.split(',')))
    
    # convert string to integer
    args.n_simulations=int(args.n_simulations)

    # Check if the directory exists and confirm with the user
    if os.path.exists(args.out_dir):
        print(f"Warning: Directory '{args.out_dir}' already exists. Proceeding may erase previous outputs.")
        if not confirm_proceed():
            exit()

    # Create output directory
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # Read tree and get taxa names
    containing_taxa = dendropy.TaxonNamespace()
    sp_tree = dendropy.Tree.get(path=args.tree, 
                                schema="nexus", 
                                preserve_underscores=True,
                                taxon_namespace=containing_taxa)

    genes_to_species = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=containing_taxa,
        num_contained=n_sampled_individuals)

    # convert to containing tree
    sp_tree = reconcile.ContainingTree(sp_tree,
                contained_taxon_namespace=genes_to_species.domain_taxon_namespace,
                contained_to_containing_taxon_map=genes_to_species)

    # Simulate and save gene trees
    trees = dendropy.TreeList()
    print('Simulating trees')
    for rep in range(args.n_simulations):
        print(rep)
        gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree, gene_to_containing_taxon_map=genes_to_species)
        trees.append(gene_tree)
    
    print('Saving newick simulated trees')
    trees.write(path= os.path.join(args.out_dir, "simulated_gene_trees.nwck"),
        schema="newick"
        )

    # Create the function to get species name of taxon object label of the gene trees. 
    # This will be used for the taxa membership, since the simulated tree has taxa with names like "species1 0", "species1 1", "species2 0"
    def mf(t):
        index=t.label.find(" ")
        return t.label[:index]
    
    # Iterate over trees to calculate 
    s_count = 0 
    target_trees = dendropy.TreeList()
    s_distribution = []
    for tree in trees:
        taxon_namespace = tree.taxon_namespace
        tax_parts = taxon_namespace.partition(membership_func=mf)
        s = monophyletic_partition_discordance(tree, taxon_namespace_partition=tax_parts)
        s_distribution.append(s)
        if s == args.observed_s_statistics:
            s_count = s_count + 1
            target_trees.append(tree)


    probability = s_count/args.n_simulations

    print('Saving newick simulated trees with target s-statistics')
    target_trees.write(path= os.path.join(args.out_dir, "target_gene_trees.nwck"),
        schema="newick"
        )

    print(f"Probability of s statistics being equal to {args.observed_s_statistics}: {probability}")

    results = open(os.path.join(args.out_dir, "S_statsProbs.txt"), "w")
    print(f"Probability of s statistics being equal to '{args.observed_s_statistics}': '{probability}'", file=results)
    results.close()

    # Save simulated s values
    df = pd.DataFrame({'simulation': range(1, len(s_distribution) + 1), 's_statistic': s_distribution})
    df.to_csv(os.path.join(args.out_dir,"simulated_s.csv"), index=False)
    print("Simulated s values saved in", os.path.join(args.out_dir,"simulated_s.csv") )

    # Create a histogram of s values
    # Compute the 95% confidence interval (2.5% and 97.5% percentiles)
    lower_bound = np.percentile(s_distribution, 5.0)

    # Create histogram
    bins = np.histogram_bin_edges(s_distribution, bins=10)
    hist_values, bin_edges, patches = plt.hist(s_distribution, bins=bins, density=True, edgecolor='black', alpha=0.7)

    # Color bars conditionally
    for patch, left_edge in zip(patches, bin_edges[:-1]):
        if left_edge < lower_bound:
            patch.set_facecolor('red')  # Outliers in red
        else:
            patch.set_facecolor('blue')  # Normal density in blue


    # Simple histogram without 95% interval
    #plt.hist(s_distribution, bins=10, edgecolor='black', alpha=0.7)

    # Add vertical line at 'target'
    plt.axvline(x=args.observed_s_statistics, color='r',
                linestyle='dashed', linewidth=2,
                label=f'Observed s = {args.observed_s_statistics}')

    # Compute and plot the smooth PDF curve with gaussian
    kde = gaussian_kde(s_distribution)  # Kernel Density Estimation
    x_vals = np.linspace(min(s_distribution), max(s_distribution), 1000)  # X range for smooth curve
    plt.plot(x_vals, kde(x_vals), color='black', linewidth=2, label="PDF Curve")  # Smooth PDF

    # Compute and plot the smooth PDF curve with lognormal distribution
    #shape, loc, scale = lognorm.fit(s_distribution, floc=0)  # Estimate parameters
    #x_vals = np.linspace(min(s_distribution), max(s_distribution), 1000)  # Smooth x values
    #pdf_vals = lognorm.pdf(x_vals, shape, loc, scale)  # Compute lognormal PDF
    ## Plot the lognormal PDF curve
    #plt.plot(x_vals, pdf_vals, color='black', linewidth=2, linestyle='dashed', label="Lognormal PDF")

    # Add Label at the Top, Close to the Line 
    y_top = plt.ylim()[1] * 1.01  # Position above the highest histogram bar
    plt.text(args.observed_s_statistics, y_top, f"p(s={args.observed_s_statistics})={probability}", color='black', ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Labels and title
    plt.xlabel('Simulated s statistics')
    plt.ylabel('Probability')

    # Save as PNG and PDF
    plt.savefig(os.path.join(args.out_dir, "histogram.png"), dpi=300)
    plt.savefig(os.path.join(args.out_dir,"histogram.pdf") )
    print("Histograms saved in", os.path.join(args.out_dir,"histogram.pdf"), "and in", os.path.join(args.out_dir,"histogram.png") )

if __name__ == "__main__":
    main()