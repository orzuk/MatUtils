load('BUGGY_EXAMPLE.mat')
use_x_flag = 1; % This flag says if we want to calculate log(Pr(X,Y)) or just log(Pr(Y))


i=22; % The last chromosome

SHORT_PLACE_M3_OFFSET = SHORT_PLACE_M3;
SHORT_PLACE_M3_OFFSET(1:end-1,:) = SHORT_PLACE_M3(2:end,:);
%SHORT_PLACE_M3_OFFSET(2:end,:) = SHORT_PLACE_M3(1:end-1,:); % Try different offset

HMM_MODEL.SIGMA2 = HMM_MODEL.SIGMA*1.1;
%HMM_MODEL.SIGMA2(:) = mean(HMM_MODEL.SIGMA2);

% Switch genotypes 0 and 1 
SHORT_PLACE_M3([30,289, 309, 346],:) = 0.5; % 346
SHORT_PLACE_M3(102,:) = [0.01 0.99 0.01 0.99];
SHORT_PLACE_M3(146,:) = [0.01 0.99 0.01 0.99];
SHORT_PLACE_M3(31,:) = [0.99 0.01 0.99 0.01];
SHORT_PLACE_M3(find(SHORT_PLACE_M3 == 0.99)) = 0.999;
SHORT_PLACE_M3(find(SHORT_PLACE_M3 < 0.011)) = 0.001;

SHORT_PLACE_M3_C = SHORT_PLACE_M3;
SHORT_PLACE_M3_C(:,1) = SHORT_PLACE_M3(:,2);
SHORT_PLACE_M3_C(:,2) = SHORT_PLACE_M3(:,1);
SHORT_PLACE_M3_C(:,3) = SHORT_PLACE_M3(:,4);
SHORT_PLACE_M3_C(:,4) = SHORT_PLACE_M3(:,3);
% Switch 2 and 3 
SHORT_PLACE_M3_C(:,2) = SHORT_PLACE_M3(:,3);
SHORT_PLACE_M3_C(:,3) = SHORT_PLACE_M3(:,2);
%TEMP = SHORT_PLACE_M3(:,2); SHORT_PLACE_M3(:,2) = SHORT_PLACE_M3(:,3); SHORT_PLACE_M3(:,3) = TEMP;

[Viterbi_Path{i} TempGammaProbs] = ...
    FindBestPathViterbi(CHROM_MATS.data_A{i}, CHROM_MATS.data_B{i}, L, ...
    use_locations, HMM_MODEL.PI, HMM_MODEL.M2, HMM_MODEL.N, ...
    HMM_MODEL.MEW', HMM_MODEL.SIGMA2, HMM_MODEL.PLACE_FLAG, ...
    HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.SNP_specific, mean_vec, std_vec, std_vec, ...
    1.0*SHORT_PLACE_M3_C'+0.0);

DataLogLikelihood = ...
    ComputeHMMLogLikelihood(Viterbi_Path{i}, CHROM_MATS.data_A{i}, CHROM_MATS.data_B{i}, use_x_flag, L, ...
    use_locations, HMM_MODEL.PI, HMM_MODEL.M2, HMM_MODEL.N, ...
    HMM_MODEL.MEW', HMM_MODEL.SIGMA2, HMM_MODEL.PLACE_FLAG, ...
    HMM_MODEL.SPECIAL_MODELS_FLAG, mean_vec, std_vec, ...
    SHORT_PLACE_M3_C')

%%                HMM_MODEL.PLACE_M3(HMM_MODEL.CHR_STARTS(i):HMM_MODEL.CHR_ENDS(i),:));  % add the place variables all_chr_exp_arr(j,:);
Viterbi_Path_joint_genotype = ( 1*mod(Viterbi_Path{i},2)+2*mod(floor(Viterbi_Path{i}/16),2) )';
HMMOutStruct.SNPsProbsMat{i} = TempGammaProbs;
[VitStruct{i} ProbsStruct{i}] = ...
    GetBestMarginalPredictions(TempGammaProbs, HMM_MODEL, do_couples, joint_flag);
sum(AlmostCorrectVitStruct.joint_genotype~=VitStruct{i}.joint_genotype)
sum(AlmostCorrectVitStruct.joint_genotype~=Viterbi_Path_joint_genotype)

sum(bitget(AlmostCorrectVitStruct.joint_genotype,1)+bitget(AlmostCorrectVitStruct.joint_genotype,2) ~= ...
    bitget(VitStruct{i}.joint_genotype,1)+bitget(VitStruct{i}.joint_genotype,2))
sum(bitget(AlmostCorrectVitStruct.joint_genotype,1)+bitget(AlmostCorrectVitStruct.joint_genotype,2) ~= ...
    bitget(Viterbi_Path_joint_genotype,1)+bitget(Viterbi_Path_joint_genotype,2))


find(bitget(AlmostCorrectVitStruct.joint_genotype,1)+bitget(AlmostCorrectVitStruct.joint_genotype,2) ~= ...
    bitget(VitStruct{i}.joint_genotype,1)+bitget(VitStruct{i}.joint_genotype,2))
f=find(bitget(AlmostCorrectVitStruct.joint_genotype,1)+bitget(AlmostCorrectVitStruct.joint_genotype,2) ~= ...
    bitget(Viterbi_Path_joint_genotype,1)+bitget(Viterbi_Path_joint_genotype,2))




ViterbiPath2 = VitStruct{i}.alpha_genotype + 16*VitStruct{i}.beta_genotype + ...
    2*VitStruct{i}.alpha_copy + 32*VitStruct{i}.beta_copy;

AlmostCorrectPath = AlmostCorrectVitStruct.alpha_genotype + 16*AlmostCorrectVitStruct.beta_genotype + ...
    2*AlmostCorrectVitStruct.alpha_copy + 32*AlmostCorrectVitStruct.beta_copy;
AlmostCorrectDataLogLikelihood = ...
    ComputeHMMLogLikelihood(AlmostCorrectPath, CHROM_MATS.data_A{i}, CHROM_MATS.data_B{i}, use_x_flag, L, ...
    use_locations, HMM_MODEL.PI, HMM_MODEL.M2, HMM_MODEL.N, ...
    HMM_MODEL.MEW', HMM_MODEL.SIGMA2, HMM_MODEL.PLACE_FLAG, ...
    HMM_MODEL.SPECIAL_MODELS_FLAG, mean_vec, std_vec, ...
    SHORT_PLACE_M3_C')



MatlabAlmostLogLikeAlpha = zeros(1,426); MatlabAlmostLogLikeBeta = zeros(1,426);
MatlabViterbiLogLikeAlpha = zeros(1,426); MatlabViterbiLogLikeBeta = zeros(1,426);
for j=1:426
    MatlabAlmostLogLikeAlpha(j) =  ...
        log(SHORT_PLACE_M3(j, 0*AlmostCorrectVitStruct.alpha_genotype(j)+1-AlmostCorrectVitStruct.alpha_genotype(j)+1));
    MatlabAlmostLogLikeBeta(j) =  ...
        log(SHORT_PLACE_M3(j, 0*AlmostCorrectVitStruct.beta_genotype(j)+1-AlmostCorrectVitStruct.beta_genotype(j)+1));
    MatlabViterbiLogLikeAlpha(j) =  ...
        log(SHORT_PLACE_M3(j, 0*mod(Viterbi_Path{i}(j),2)   + 1-mod(Viterbi_Path{i}(j),2) +1));
    MatlabViterbiLogLikeBeta(j) =  ...
        log(SHORT_PLACE_M3(j, 0*mod( floor(Viterbi_Path{i}(j)/16) ,2)   + 1-mod(  floor(Viterbi_Path{i}(j)/16) ,2) +1));

end
MatlabAlmostLogLike = MatlabAlmostLogLikeAlpha  + MatlabAlmostLogLikeBeta; 
MatlabViterbiLogLike = MatlabViterbiLogLikeAlpha + MatlabViterbiLogLikeBeta;

Almost10 = MatlabAlmostLogLike(1:10)
Vit10 = MatlabViterbiLogLike(1:10)

R = [];
[R.type R.chr R.start R.stop R.ind] = StratifyToGenomicRegions([2 3], [1000000 10000], [2000000 999999], 'hg18'); % get regions

R = sort_struct_by_field(R, 'start');

intervals_plot(R.start, R.stop, R.type, 'b');





