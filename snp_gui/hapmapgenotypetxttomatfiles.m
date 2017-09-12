% Written by Or Zuk 5/2007
%
% This function transfers the .txt files containing the HAPMAP genotype data to .mat files.
% There are two possibilities : One .mat file for each chromosome OR One
% huge .mat file for the whole data. 
% Note: Currently SaveOneFile is off - this means that we don't save one files for
% all the chromosomes but save seperate files for each chromosome.
%
% The input: 
% population - which hapmap population to take
% hapmap_version - current version of the hapmap project
% chroms - which chromosomes to run on
% save_flag - 1 save into file. 0 - don't save
%
% The output:
%
% SnpsData - The data matrix itself
function [SnpsData, SnpsBases, SnpsFreqs, SnpsHetros, SnpsBadCalls, SnpsChromLocs SnpsNames] = ...
    HapMapGenotypeTxtToMatFiles(population, hapmap_version, hapmap_dir, chroms, file_names, save_flag)

AssignAllGlobalConstants();

cur_dir = pwd; cd(fullfile(hapmap_dir, pop_str_vec{population}, hapmap_version)); SaveOneFile=0;

for chrom=chroms
    do_chrom = chrom
    chr_str = chr_num2str(chrom);

%%    cur_file_name = ['genotypes_chr' num2str(chrom) '_' pop_str_vec{population} '_' hapmap_version '_nr_fwd.mat'];     
    for j=1:length(file_names)
        if(~isempty(findstr(['_chr' chr_str '_'], file_names{j})) );
            cur_file_name = file_names{j};
        end
    end

    % SnpsData is binary packed and holds the SNPs values for each individual. 0 is the major allele and 1 the minor allele
    % So we (should) have: AA 0 AB 1 BB 2 NN 3 
    % SnpsBases holds which bases are we talking about: 4 A/T, 6 C/T, 7 G/T, 8 A/C, 12 A/G, 14 C/G
    if(SaveOneFile)
        [SnpsData{chrom}, SnpsBases{chrom}, SnpsFreqs{chrom}, SnpsHetros{chrom}, ...
            SnpsBadCalls{chrom}, SnpsChromLocs{chrom} SnpsNames{chrom}] = ...
            ReadGenotypeFile(cur_file_name); %% chrom, pop_str_vec{population}, hapmap_version);
        SnpsData{chrom} = SnpsData{chrom}'; SnpsBadCalls{chrom} = SnpsBadCalls{chrom}';
    else
        [SnpsData, SnpsBases, SnpsFreqs, SnpsHetros, SnpsBadCalls, SnpsChromLocs SnpsNames] = ...
            ReadGenotypeFile(cur_file_name); %% chrom, pop_str_vec{population}, hapmap_version);
        SnpsData = SnpsData'; SnpsBadCalls = SnpsBadCalls';
        if(save_flag)
            save([cur_file_name(1:end-3) 'mat'], 'SnpsData', 'SnpsBases', 'SnpsFreqs', 'SnpsHetros', 'SnpsBadCalls', 'SnpsChromLocs', 'SnpsNames');
%%            eval(['save genotypes_chr' ...
%%                chr_str '_' pop_str_vec{population} '_' hapmap_version '_nr_fwd.mat SnpsData SnpsBases SnpsFreqs SnpsHetros SnpsBadCalls SnpsChromLocs SnpsNames;']);
        end
    end

    % Save to be sure ...
    %save ['./Database/LD_' chip_type] LD

end

if(save_flag)
    if(SaveOneFile) % Here WE determin ethe file name - not based on what is in Hapmap
        save(['All_genotypes_chr_CEU_' hapmap_version '_nr_fwd.mat'], ...
        'SnpsData', 'SnpsBases', 'SnpsFreqs', 'SnpsHetros', 'SnpsBadCalls', 'SnpsChromLocs', 'SnpsNames');
    end
end    

eval(['cd ' cur_dir]); % Go back to src directory

Dummy=99;
