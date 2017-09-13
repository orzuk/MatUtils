Genetic Architecture 
====================

This package is a sub-package of the MatUtils general package, 
and contains functions for simulation and analysis of genetic architectures of complex traits. 
It acompanies the papers: 

[1] "The Mystery of Missing Heritability: Genetic Interactions Create Phantom Heritability" 
O. Zuk, E. Hechter, S. Sunyaev and E.S. Lander 
Proceedings of the National Academy of Sciences, 109(4): 1193-1198  (2012) 

[2] "Searching for Missing Heritability:  Designing Rare Variants Association Studies"
O. Zuk, S. Schaffner, K. Samocha, R. Do, E. Hechter, S. Katherisan, M. Daly, B. Neale, S.R. Sunyaev and E.S. Lander
Proceedings of the National Academy of Sciences, doi: 10.1073/pnas.1322563111 (2014) 

Please cite these papers if you're using the package. 

The package replaces and extends upon the old 'heritability calculator' (hc) package available at http://software.broadinstitute.org/mpg/hc/. 
Please don't use the hc package as it is obsolete.


Details: 
--------

The package enables heritability calculations using different mehotds and for many differnet disease models (genetic architectures)

It enables simulations of different sequencing pooling experiments, including pooling design, modeling sequencing errors, reconstruction of genotype vectors and evaluating reconstruction accuracy. 
It may serve as a guideline for designing pooling sequencing experiments.
One may also enter real allele counts collected in pooled experiments and reconstruct the genotype.

Please look at the following files for the different applications:

(1) **test_estimate_heritability_using_IBD.m** - A script testing heritability estimation by local-regression on IBD. The script simulates IBD structure, genotypes 
and phenotypes, with phenotypes being a non-linear noisy function of genotypes. Then, we estimate narrow-sense heritability by regressing phenotypic correlation 
on IBD sharing and compare to the true narrow-sense heritability.

(2) compute_k_of_N_liabilities_statistics.m - computing all relevant genetic epidemiological parameters for a generalized Limiting-Pathway genetic architecture.

(3) test_estimate_heritability_using_IBD.m - Simulate genotypes for a population with common founders (shared IBD blocks), simulate phenotypes,
 and estimated heritability by estimating of IBD sharing between pairs of individuals (used for heritability calculations).

(4) compute_association_power.m - Power calculations for different tests (e.g. marginal effect/epistasis/pathway-epistasis tests) and different genetic architectures. 

(5) RVAS_power_calculator.m - Power calculation for aggregate tests for rare variants. Can also compute power as function of population genetics parameters (e.g. selection)

(6) plot_two_class_equilibrium_statistics.m, plot_two_class_empirical_distributions.m - compute and plot several statistics of the allelic frequency spectrum as function of population genetics parameters 

(7) plot_RVAS_paper_figs.m - Plot several of the figures appearing in [2]. 


(*) The package was not used to perform forward simulations for populations which appeared in [2] to get simulated allelic data. Such simulations were performed by Steven Schaffner (sfs@broadinstitute.org)

(*) The package was tested on Matlab 2016a, but should work on most Matlab versions. 

Please send any comments/bug reports to: Or Zuk: or.zuk@mail.huji.ac.il


Ver 4.0: Jan, 2017

