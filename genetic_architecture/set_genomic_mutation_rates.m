% Set parameters used for mutation rates for typical genes 
%
% Input:
% alpha - fraction of nulls at birth 
% beta - fraction of hypermorphs at birth 
% 
% Output:
% mu - total mutation rate per-nucleotide 
% gene_length - number of codons in total exons of typical gene
% mu_pwer_gene - strucutre with mutation rates for different types of mutations
% 
function [mu gene_length mu_per_gene] = set_genomic_mutation_rates(alpha, beta)

if(~exist('alpha', 'var') || isempty(alpha))
    alpha = 0.25; % take default
end
if(~exist('beta', 'var') || isempty(beta))
    beta = 0; % take default: no hypermorphs
end

mu=3.077*10^(-8); % why so high?? 1.6*10^(-8); % mutation rate per-nucleotide per-generation
gene_length = 625; % 1500; % typical gene length in nucleotides. Makes the gene mutation rate mu_g at 10^(-5)
mu_per_gene.NULL = mu*gene_length; % Set some parameters for different classes
mu_per_gene.LOF = mu*gene_length * 0.047; mu_per_gene.NONSENSE = mu_per_gene.LOF;
mu_per_gene.FRAMESHIFT = mu*gene_length * 0.043;
mu_per_gene.DISRUPTIVE = mu_per_gene.LOF + mu_per_gene.FRAMESHIFT;
mu_per_gene.MISSENSE = mu*gene_length * 0.663; % (4/5)*2; % assume alpha=1/2
mu_per_gene.MISSENSE_NULL = mu_per_gene.MISSENSE * alpha;
mu_per_gene.MISSENSE_HYPERMORPH = mu_per_gene.MISSENSE * beta;
mu_per_gene.MISSENSE_SILENT = mu_per_gene.MISSENSE * (1-alpha-beta); mu_per_gene.MISSENSE_NEUTRAL = mu_per_gene.MISSENSE_SILENT;

mu_per_gene.SYNONYMOUS = mu*gene_length * 0.29;
mu_per_gene.TOTAL_SUBSTITUTIONS = mu_per_gene.LOF+mu_per_gene.MISSENSE+mu_per_gene.SYNONYMOUS; % total mutation rate (should include frame-shifts)
mu_per_gene.TOTAL = mu_per_gene.TOTAL_SUBSTITUTIONS + mu_per_gene.FRAMESHIFT; 
