% Written by Or Zuk 5/2007
%
% This function reads the genotype data of ALL samples in the HAPMAP project, 
% and takes only the SNPs which appear on a given chip, and saves this data to another file
% in the database
function Dummy = HapMapExtractOnlyChipGenotypes(chroms, chip_type, hapmap_version, hapmap_dir, hapmap_population, file_names)


% load the SNPs which appear on the chip (e.g. hind, xba or both) ...
AssignAllGlobalConstants();

% Old: Put CHIPs together: 
% load(['../../DataBase/Hind_annot_data_' genome_assembly '.mat'], 'snp_ids', 'rs_ids', 'chr_vec'); snp_hind_ids = snp_ids; rs_hind_ids = rs_ids; chr_hind_vec = chr_vec; 
% load(['../../DataBase/Xba_annot_data_' genome_assembly '.mat'], 'snp_ids', 'rs_ids', 'chr_vec'); snp_xba_ids = snp_ids; rs_xba_ids = rs_ids; chr_xba_vec = chr_vec;
% snp_ids = [snp_hind_ids' snp_xba_ids']'; rs_ids = [rs_hind_ids' rs_xba_ids']'; chr_vec = [chr_hind_vec' chr_xba_vec']';


% One file for each CHIP
load(fullfile('../DataBase', [chip_type '_annot_data_' genome_assembly '.mat']), 'snp_ids', 'rs_ids', 'chr_vec'); 
rs_ids_char = char(rs_ids); 

cur_dir = pwd; cd(fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version)); %   E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\CEU; 


% Go chrom-by-chrom as to not take too much memory
for chrom=chroms
    do_chrom = chrom
    chr_str = chr_num2str(chrom); % get the chromosome string
    
    
    % load the file 
    for j=1:length(file_names)
        if(~isempty(findstr(['_chr' chr_str '_'], file_names{j})) );
            cur_file_name = file_names{j};
        end
    end    
    
    R = load(cur_file_name);    
 %%   eval(['R=load(''genotypes_chr' chr_str '_' pop_str_vec{hapmap_population} '_' hapmap_version '_nr_fwd.mat'');']);
    R.SnpsNames(R.SnpsNames == char(0)) = ' '; % Convert 0 chars (null?) to spaces
    
    [rs_intersect I J] = intersect(rs_ids_char, R.SnpsNames, 'rows');
    
    SnpsNames{chrom} = R.SnpsNames(J,:); % Extract only the common SNPs !! 
    SnpsData{chrom} = R.SnpsData(J,:);
    SnpsBases{chrom} = R.SnpsBases(J);
    SnpsFreqs{chrom} = R.SnpsFreqs(J);
    SnpsHetros{chrom} = R.SnpsHetros(J);
    SnpsBadCalls{chrom} = R.SnpsBadCalls(J,:);
    SnpsChromLocs{chrom} = R.SnpsChromLocs(J);
    SnpsIDs{chrom} = snp_ids(I); % Here save the SNPs (not the rs) IDS

end

SaveOneFile=1;
if(SaveOneFile)
    save([chip_type '_genotypes_chr_' pop_str_vec{hapmap_population} '_' hapmap_version '_nr_fwd.mat'], ...
        'SnpsData', 'SnpsBases', 'SnpsFreqs', 'SnpsHetros', 'SnpsBadCalls', 'SnpsChromLocs', 'SnpsNames', 'SnpsIDs');
%%    HindXba_genotypes_chr_CEU_r21_nr_fwd.mat, 'SnpsData', 'SnpsBases', 'SnpsFreqs', 'SnpsHetros', ...

end

eval(['cd ' cur_dir]); % Go back to src directory

Dummy=99;
