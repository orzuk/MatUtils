% Create and input according to the format of the Displaying function we use
function [SampleChromsDataStruct] = PrepareInputForChromsDislayFunction(Gamma_Probs, HMM_MODEL_SAVED, HMM_chrom_loc_mat)

do_couples = 1; HMM_x_dim = 3;
Upperbound = 2.5; Lowerbound = 1.5;

SampleChromsDataStruct = {};

% first create empty matrices
for i=1:24
    SampleChromsDataStruct{i} = [];
end


% Compute data for all otozomal chromosomes
for i=1:22
    I_is = i
    %     Gamma_Probs
    %      HMM_MODEL_SAVED{i}{HMM_x_dim}
% %     [VitStruct.alpha_genotype VitStruct.beta_genotype VitStruct.alpha_copy VitStruct.beta_copy ...
% %         VitStruct.total_copy VitStruct.A_copy VitStruct.B_copy] = ...
     [VitStruct ProbsStruct] = ...    
        GetBestMarginalPredictions(reshape(Gamma_Probs{i}(1,:,:), 256, size(Gamma_Probs{i}, 3) ), ...
        HMM_MODEL_SAVED{i}{HMM_x_dim}, do_couples);

    SampleChromsDataStruct{i}.Data = [VitStruct.alpha_copy VitStruct.beta_copy];
    SampleChromsDataStruct{i}.Locs = HMM_chrom_loc_mat{i};


    AmpInds = (VitStruct.alpha_copy+VitStruct.beta_copy > Upperbound)'; % Find the segments positions
    [AmpIndsStarts AmpIndsEnds] = GetSegmentsFromIndexes(AmpInds);
    DelInds = (VitStruct.alpha_copy+VitStruct.beta_copy < Lowerbound)';
    [DelIndsStarts DelIndsEnds] = GetSegmentsFromIndexes(DelInds);
    AllIndsStarts = [AmpIndsStarts DelIndsStarts];
    AllIndsEnds = [AmpIndsEnds DelIndsEnds];
    [AllIndsStarts SortPerm] = sort(AllIndsStarts); AllIndsEnds = AllIndsEnds(SortPerm); % Sort segments

    SampleChromsDataStruct{i}.Segments(:,1) = HMM_chrom_loc_mat{i}(AllIndsStarts);
    SampleChromsDataStruct{i}.Segments(:,2) = HMM_chrom_loc_mat{i}(AllIndsEnds);
    V_VEC_cumsum = [0 cumsum(VitStruct.alpha_copy+VitStruct.beta_copy)'];
    SampleChromsDataStruct{i}.Segments(:,3) = (V_VEC_cumsum(AllIndsEnds+1) - V_VEC_cumsum(AllIndsStarts)) ./ ...
        (AllIndsEnds-AllIndsStarts+1); % Find the average of each semgment
end







