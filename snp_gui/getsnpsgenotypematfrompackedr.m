function GenotypeMat = GetSNPsGenotypeMatFromPackedR(R, snp_ids); 
AssignAllGlobalConstants;
chroms = 1:24;

num_samples = 90; 

GenotypeMat = zeros(num_samples, length(snp_ids));

SampleInds = 1:num_samples; % samples order doesn't matter 
SampleBits = mod(SampleInds-1,30)+1; SampleWords = ceil(SampleInds./30);


for chrom=chroms
    chrom_is = chrom
    [snp_ids_intersect I J] = intersect(snp_ids, R.SnpsIDs{chrom}); % Get the common SNPs from disp

    CurSNPsDataChromMat = R.SnpsData{chrom}(J,:);
    CurSNPsBadCallChrom = find(bitget(R.SnpsBadCalls{chrom}(J,1),1));

    for cur_sample = 1:num_samples
        sample_is = cur_sample

        HAPMAP_SNPs = 2*bitget(CurSNPsDataChromMat(:,SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
            bitget(CurSNPsDataChromMat(:,SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP. Bad Calls are considered as '-1'
        HAPMAP_SNPs(CurSNPsBadCallChrom) = -1; % Bad calls NN SNPs are set to minus one

        HAPMAP_SNPs(HAPMAP_SNPs == 2) = 1; % unite AB and BA

        % Transfer from (0,1,3,-1) to (1,2,3,4)
        HAPMAP_SNPs(HAPMAP_SNPs == 3) = BB;
        HAPMAP_SNPs(HAPMAP_SNPs == 1) = AB;
        HAPMAP_SNPs(HAPMAP_SNPs == 0) = AA;
        HAPMAP_SNPs(HAPMAP_SNPs == -1) = NoCall;

        GenotypeMat(cur_sample,I) = uint8(HAPMAP_SNPs);
    end
end

GenotypeMat = GenotypeMat';

