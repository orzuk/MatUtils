% This script is used for updating the LD structure 
use_hetro_flag = 0; % We currently don't know how to use them
chroms = 1:24; % all otozomal chromosomes, and also X and Y
%chip_type = 'Xba';
AssignAllGlobalConstants(); 
hapmap_population = CEU; hapmap_version = 'r21';
hapmap_dir = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\';

LD_Xba = ComputeSNPsAdjacentLD(chroms, hapmap_dir, 'Xba', hapmap_population, hapmap_version, genome_assembly, use_hetro_flag);
LD_Hind = ComputeSNPsAdjacentLD(chroms, hapmap_dir, 'Hind', hapmap_population, hapmap_version, genome_assembly, use_hetro_flag);