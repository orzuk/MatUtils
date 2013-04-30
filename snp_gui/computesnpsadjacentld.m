% Written by Or Zuk 6/2007
%
% This function computes the Linkage Disequilibrium values for pairs of SNPs which are on the
% Affy chips. The function also saves the resulting LD structure into files.
%
% The inputs are:
% chroms - which chromsomes to do. Input should be all chromosomes (1-24).
% chip_type - type of chip (e.g. 'xba')
% genome_assembly - genome version (e.h. 'hg18')
% use_heter_flag - flag saying if to use heterozygous SNPs for statistics (currently not used)
%
% Output: 
% LD - a structure with Linkage disequilibrium statistics
function LD = ComputeSNPsAdjacentLD(chroms, hapmap_dir, chip_type, hapmap_population, hapmap_version, genome_assembly, use_hetro_flag)
AssignAllGlobalConstants();

% Load chip files
ChipStr = [chip_type '_annot_data_' genome_assembly '.mat'];
load(fullfile('..', 'Database', ChipStr), 'rs_ids', 'chr_loc_vec');


Calc_LD_mat_flag = 0;

% first read the xba chip file
%SamplesInfo = loadCellfile(fullfile(user_dir, 'E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\pedinfo2sample_CEU.txt');
SamplesInfo = loadCellfile(fullfile(hapmap_dir, pop_str_vec{hapmap_population} , hapmap_version, ...
    ['pedinfo2sample_' pop_str_vec{hapmap_population} '.txt']));

SamplesFamilyID = str2num(cell2mat(SamplesInfo(:,1)));
SamplesLabels =  SamplesInfo(:,7);
for i=1:90
    SamplesLabels{i} = SamplesLabels{i}(end-8:end-2);
end
ParentsLabels = SamplesLabels(setdiff(1:90, 3:3:90)); OffspringsLabels = SamplesLabels(3:3:90);

% The same samples, different order ..
[SamplesFromChromLabels SortInds] = sort(SamplesLabels);
[OffspringsIntersect OffspringsInds OffspringsJnds] = intersect(SamplesFromChromLabels, OffspringsLabels);
%% WronGGGGGGGG !!!!  OffspringsInds = SortInds(3:3:90);

%% NA06985 NA06991 NA06993 NA06994 NA07000 NA07019 NA07022 NA07029 NA07034 NA07048 NA07055 NA07056 NA07345 NA07348 NA07357 NA10830 NA10831 NA10835 NA10838 NA10839 NA10846 NA10847 NA10851 NA10854 NA10855 NA10856 NA10857 NA10859 NA10860 NA10861 NA10863 NA11829 NA11830 NA11831 NA11832 NA11839 NA11840 NA11881 NA11882 NA11992 NA11993 NA11994 NA11995 NA12003 NA12004 NA12005 NA12006 NA12043 NA12044 NA12056 NA12057 NA12144 NA12145 NA12146 NA12154 NA12155 NA12156 NA12234 NA12236 NA12239 NA12248 NA12249 NA12264 NA12707 NA12716 NA12717 NA12740 NA12750 NA12751 NA12752 NA12753 NA12760 NA12761 NA12762 NA12763 NA12801 NA12802 NA12812 NA12813 NA12814 NA12815 NA12864 NA12865 NA12872 NA12873 NA12874 NA12875 NA12878 NA12891 NA12892

%% Big File Ordering:
%%%%%%%%%%%%%%%%%%%%%
%% NA06985 NA06991 NA06993 NA06994 NA07000 NA07019 NA07022 NA07029 NA07034
%% NA07048 NA07055 NA07056 NA07345 NA07348 NA07357 NA10830 NA10831 NA10835
%% NA10838 NA10839 NA10846 NA10847 NA10851 NA10854 NA10855 NA10856 NA10857
%% NA10859 NA10860 NA10861 NA10863 NA11829 NA11830 NA11831 NA11832 NA11839
%% NA11840 NA11881 NA11882 NA11992 NA11993 NA11994 NA11995 NA12003 NA12004
%% NA12005 NA12006 NA12043 NA12044 NA12056 NA12057 NA12144 NA12145 NA12146
%% NA12154 NA12155 NA12156 NA12234 NA12236 NA12239 NA12248 NA12249 NA12264
%% NA12707 NA12716 NA12717 NA12740 NA12750 NA12751 NA12752 NA12753 NA12760
%% NA12761 NA12762 NA12763 NA12801 NA12802 NA12812 NA12813 NA12814 NA12815
%% NA12864 NA12865 NA12872 NA12873 NA12874 NA12875 NA12878 NA12891 NA12892




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start LD computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The linkage disequilibrium struct.
% Can contain several different measures

% Use C function to load the hapmap data for all the chromosomes (which C function ??)
cur_dir = pwd; path(path, cur_dir);
for chrom=chroms
    do_chrom = chrom
    
    cur_dir = pwd; cd(fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version));
%%    cd (fullfile(hapmap_dir, 'E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\CEU;
    % SnpsData is binary packed and holds the SNPs values for each individual. 1 is the major allele and 0 the minor allele
    % SnpsBases holds which bases are we talking about: 4 A/T, 6 C/T, 7 G/T, 8 A/C, 12 A/G, 14 C/G    
%%    cur_file_name = ['genotypes_chr' num2str(chrom) '_' 'CEU' '_' 'r21' '_nr_fwd.mat']; 
%%    cur_file_name = ['genotypes_chr' chr_num2str(chrom) '_' pop_str_vec{hapmap_population} '_' hapmap_version '_nr_fwd.mat']; 

    s = dir(['*chr' chr_num2str(chrom) '_*.mat']); cur_file_name = s.name; % get file name from directory
%%    [SnpsData, SnpsBases, SnpsFreqs, SnpsHetros, SnpsBadCalls, SnpsChromLocs SnpsNames] = ReadGenotypeFile(cur_file_name);
    load(cur_file_name); % since it's already in .mat format we just need to load it ..
    eval(['cd ' cur_dir]);
    SnpsData = SnpsData'; SnpsBadCalls = SnpsBadCalls';

    % Find which SNPs are on the input chip AND on the chromosome
    [names I J] = (intersect(rs_ids, SnpsNames(2:end,:)));  J=J+1; % I=I+1;
    [Vals Inds] = sort(chr_loc_vec(I));
    names = names(Inds); I = I(Inds); J = J(Inds);
    % Call the function to calculate all the correlations. We ignore the
    % offsprings when doing the statistics!!!  If we do want to insert
    % offsprings, give an empty vector for the OffspringsInds
    LD{chrom}.RS_IDs = names; 
    % If we want the entrire LD pairwise matrix we should replace r_mat by LD{chrom}.r_mat
    [r_mat LD{chrom}.FreqVec LD{chrom}.PairMat]= ...
        CalcLinkageDisequilirium(SnpsFreqs(J), SnpsData(:,J)', SnpsBadCalls(:,J)', OffspringsInds, 1, 1, Calc_LD_mat_flag); % Don't use Hetros

    save(fullfile('..', 'Database', ['LD_' pop_str_vec{hapmap_population} '_' chip_type '_' genome_assembly]), 'LD');     % Save to be sure  at each chromsome ... 
end

