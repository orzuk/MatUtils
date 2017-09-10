% Compute Linkage Disequilibrium values for pairs of SNPs which are on the
% Affy chips!!!
function LD = ComputeSNPsAdjacentLD(chroms, use_hetro_flag)
% Don't work twice  ... 
load('LD_all_chroms.mat');  %%LD = {};

Calc_LD_mat_flag = 0;
X=200; % Number of first SNPs to take

% first read the xba chip file
SamplesInfo = loadcellfile('E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\pedinfo2sample_CEU.txt');
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


% Load chip files
load('Xba_allele_info_table.mat');
did_Xba = 999
load('Hind_allele_info_table.mat');
did_Hind = 888


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start LD computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The linkage disequilibrium struct.
% Can contain several different measures

% Use C function to load the hapmap data for all the chromosomes
for chrom=chroms
    do_chrom = chrom
    
    [SnpsData, SnpsBases, SnpsFreqs, SnpsHetros, SnpsBadCalls, SnpsChromLocs SnpsNames] = ReadGenotypeFile(chrom);
    SnpsData = SnpsData';SnpsBadCalls = SnpsBadCalls';

    % Find which SNPs are on the Xba chip AND on the chromosome
    [names I J] = (intersect(xba_table(2:end,4), SnpsNames(2:end,:))); I=I+1; J=J+1;
    [Vals Inds] = sort(str2num(char(xba_table(I,3))));
    names = names(Inds); I = I(Inds); J = J(Inds);
    % Call the function to calculate all the correlations. We ignore the
    % offsprings when doing the statistics!!!  If we do want to insert
    % offsprings, give an empty vector for the OffspringsInds
    LD{chrom}.Xba.Names = names; LD{chrom}.Xba.ChromLocs = SnpsChromLocs(J);
    [LD{chrom}.Xba.r_mat LD{chrom}.Xba.FreqVec LD{chrom}.Xba.PairMat]= ...
        CalcLinkageDisequilirium(SnpsFreqs(J), SnpsData(J,:), SnpsBadCalls(J,:), OffspringsInds, 1, 1, Calc_LD_mat_flag); % Don't use Hetros

    % Find which SNPs are on the Hind chip AND on the chromosome
    [names I J] = (intersect(hind_table(2:end,4), SnpsNames(2:end,:))); I=I+1; J=J+1;
    [Vals Inds] = sort(str2num(char(hind_table(I,3))));
    names = names(Inds); I = I(Inds); J = J(Inds);
    % Call the function to calculate all the correlations. We ignore the
    % offsprings when doing the statistics!!!  If we do want to insert
    % offsprings, give an empty vector for the OffspringsInds
    LD{chrom}.Hind.Names = names; LD{chrom}.Hind.ChromLocs = SnpsChromLocs(J);
    [LD{chrom}.Hind.r_mat LD{chrom}.Hind.FreqVec LD{chrom}.Hind.PairMat]= ...
        CalcLinkageDisequilirium(SnpsFreqs(J), SnpsData(J,:), SnpsBadCalls(J,:), OffspringsInds, 1, 1, Calc_LD_mat_flag); % Don't use Hetros    

    % Save to be sure ... 
    save 'LD_all_chroms' LD

end

save 'LD_all_chroms' LD

% Compare SNPs frequencies
%figure; plot( min(str2num(cell2mat(xba_table(I,9))),1-str2num(cell2mat(xba_table(I,9)))), min(SnpsFreqs(J),1-SnpsFreqs(J)), '.');
%xlabel('Affy HAPMAP Freqs'); ylabel('My Freqs from C');

%%[LD.r_mat1 LD.FreqVec]= CalcLinkageDisequilirium(SnpsFreqs(J), SnpsData(J,:), SnpsBadCalls(J,:), OffspringsInds, 1, 1); % Don't Use Hetros
%figure; plot(LD.r_mat1(:), LD.r_mat(:), '.'); title('My comparison'); xlabel('Use Hetros'); ylabel('Dont use Hetros');
%figure; plot( min(str2num(cell2mat(xba_table(I,9))),1-str2num(cell2mat(xba_table(I,9)))), min(LD.FreqVec,1-LD.FreqVec), '.');
%xlabel('Affy HAPMAP Freqs'); ylabel('My Freqs from Matlab');
%figure; imagesc(LD.r_mat.^2); colorbar; title('Xba SNPs LD');


% % % % % % Just for checking .. This part is removed for now
% % % % % % Here take consecutive hapmap SNPs (not only from the CHIP)
% % % % % [LD.FirstXSNPsMat LD.FirstXLocalFreqVec BC0]= ...
% % % % %     CalcLinkageDisequilirium(SnpsFreqs(2:X+1), SnpsData(2:X+1,:), SnpsBadCalls(2:X+1,:), OffspringsInds, 1, 1);
% % % % % [LD.FirstXSNPsMat1 LD.FirstXLocalFreqVec1 BC1]= ...
% % % % %     CalcLinkageDisequilirium(SnpsFreqs(2:X+1), SnpsData(2:X+1,:), SnpsBadCalls(2:X+1,:), OffspringsInds, 1, 1);
% % % % % figure; plot(LD.FirstXSNPsMat1(:), LD.FirstXSNPsMat(:), '.'); title('My comparison'); xlabel('Use Hetros r^2'); ylabel('Dont use Hetros r^2');
% % % % %
% % % % %
% % % % % figure; imagesc(LD.FirstXSNPsMat.^2); colorbar; title('Consecutive first ??? SNPs');
% % % % %
% % % % % % Check if we and HAPMAP get the same LD r^2 values:
% % % % % HAPMAP_firstXFreqs = loadcellfile('E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\CHR22_First200_Freqs.txt');
% % % % % HAPMAP_freq_vec = str2num(char(HAPMAP_firstXFreqs(4:end,7)));
% % % % % % Get intersection
% % % % % [VV III JJJ]  = intersect(HAPMAP_firstXFreqs(4:end,1), SnpsNames(2:X+1,:)); % First SNP
% % % % % figure; plot(min(HAPMAP_freq_vec,1-HAPMAP_freq_vec), min(LD.FirstXLocalFreqVec(JJJ),1-LD.FirstXLocalFreqVec(JJJ)), '.');
% % % % % xlabel('HAPMAP Freqs'); ylabel('MY Freqs');
% % % % %
% % % % %
% % % % %
% % % % % HAPMAP_firstX = loadcellfile('E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\CHR22_First200x.txt');
% % % % % HAPMAP_r_sqr_vec = str2num(char(HAPMAP_firstX(3:end-2,7)));
% % % % %
% % % % % my_r_sqr_vec = zeros(1, length(HAPMAP_r_sqr_vec));
% % % % %
% % % % % % Loop over all couples which appear in the HAPMAP vec
% % % % % for ind=3:(length(HAPMAP_firstX)-2)
% % % % %    [V I1 I2]  = intersect(HAPMAP_firstX(ind,4), SnpsNames(2:X+1,:)); % First SNP
% % % % %    [V J1 J2]  = intersect(HAPMAP_firstX(ind,5), SnpsNames(2:X+1,:)); % Second SNP
% % % % %    my_r_sqr_vec(ind-2) =  LD.FirstXSNPsMat(I2,J2)^2;
% % % % % end
% % % % %
% % % % % figure; plot(HAPMAP_r_sqr_vec, my_r_sqr_vec, '.'); xlabel('HAPMAP r^2'); ylabel('MY r^2');
% % % % %
% % % % % equal_pairs = find( (HAPMAP_r_sqr_vec - my_r_sqr_vec').^2 < 0.000001)
% % % % % figure; plot(HAPMAP_r_sqr_vec(equal_pairs), my_r_sqr_vec(equal_pairs), '.'); xlabel('HAPMAP r^2'); ylabel('MY r^2'); title('only equal');
% % % % %
% % % % %
% % % % %
