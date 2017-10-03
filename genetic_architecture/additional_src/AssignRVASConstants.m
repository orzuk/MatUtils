% Set constants and paths for RVAS project
exome_data = 'ExAC'; % 'ESP'; % 'ExAC'; % NEW! add also Exome Aggregation Data!!!!

if(~exist('machine', 'var'))
    AssignGeneralConstants;
end
if(machine == UNIX) % set sfs_code path
    % sfs_code_dir = 'c/Users/user/Downloads/sfscode';
else % PC
    sfs_code_dir = ['C:/Users/' user_str '/Downloads/sfscode'];
end
sfs_figs_dir = ['C:\Users\' user_str '\Dropbox\rare_alleles_paper\JamesZou\figs'];

exome_data_figs_dir = ['C:\Users\' user_str '\Dropbox\rare_alleles_paper\RVAS_GeneSpecific\figs\' exome_data];



switch machine % Get directory of data 
    case UNIX        
        spectrum_data_dir = '/seq/orzuk/common_disease_model/data/SiteFrequencySpectra/'; % broad
        if(~exist(specturm_data_dir, 'dir'))
            spectrum_data_dir = '/cs/cbio/orzuk/projects/RVAS/data/SiteFrequencySpectrum'; % huji
        end
        if(~exist('in_matlab_flag', 'var'))
            in_matlab_flag = 0;
        end
        if(~exist('mem_flag', 'var')) % || isempty(mem_flag))
            mem_flag = 8; % allow large memory
        end
    case PC
        %        spectrum_data_dir = 'C:\\research\common_disease_model\data\SiteFrequencySpectra'; % OLD PC
        %%%        spectrum_data_dir = 'D:\research\RVAS\Data\SiteFrequencySpectra'; % NEW PC DELL
        spectrum_data_dir = 'C:\research\RVAS\Data\SiteFrequencySpectra'; % NEW SURFACE
        
        %        spectrum_data_dir = 'T:\\common_disease_model\data\SiteFrequencySpectra';
        in_matlab_flag = 1;
end
exons_file = 'hg18_exons.mat'; % where to save exons
triplet_mutations_file = 'scone_hg17_m3_64x64_rates.txt'; % 'triplets_file_human_chimp_baboon.mat'; % file with mutation rates from all 64 codons (estimated from human-chimp-baboon alignment)
mutation_rates_file = 'human_genes_mutation_rates_hg18.mat'; % output file with triplet mutation rates

MISSENSE_C = 1;  NEUTRAL_C = 2; SYNONYMOUS_C = 2;  NULL_C = 3; % set allele classes. Missense BEFORE(!) Synonymous !! 

N_eff = 10000; % set effective population size 
mu_per_site = 2*10^(-8); % set mutation rate per-generation per nucleotide 

global cumsum_log_vec;
cumsum_log_vec = cumsum([0 log(1:2*N_eff)]);

