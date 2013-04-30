hapmap_version = genome_assembly_to_hapmap_version(); 
hapmap_dir = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data'; 
chip_types{1} = 'hind'; chip_types{2} = 'xba'; 
%hapmap_data_http = 'http://www.hapmap.org/genotypes/latest_ncbi_build36/fwd_strand/non-redundant/';

UpdateDatabaseHapmap(hapmap_version, hapmap_dir, chip_types);


