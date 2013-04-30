% Written by Or Zuk 6/2007
%
% This function transfers the linkage disequilibrium statistics into
% place-dependent Markov transition probability matrices.
% In addition, we save also the locations of chromosomes starts and ends.
%
% The Inputs are: 
% LDStruct - structure containing Linkeage-Disequilibrium statistics
% HMM_MODEL_IN - structure containing input model parameters
% HMMParamsStruct - structure contining input parameters
% 
% Output: 
% HMM_MODEL - structure containing output model parameters
%
function HMM_MODEL = LDToPlaceM(LDStruct, HMM_MODEL_IN, HMMParamsStruct)

HMM_MODEL = HMM_MODEL_IN;

% First loop just to determine size
cur_ind=1;
for i=HMMParamsStruct.ChromosomesToRun
    cur_ind = cur_ind+length(LDStruct.LD{i}.PairMat(:,1))+1; % Take extra one, skip ..
end

% Always update the PLACE_M according to the whole data
HMM_MODEL.PLACE_M = zeros(cur_ind-1,4); HMM_MODEL.CHR_STARTS = zeros(1,24); HMM_MODEL.CHR_ENDS = zeros(1,24); 
cur_ind=1;

switch HMMParamsStruct.UseGenotypesCorrs 
    case 0 % use 'full' LD pairwise statsitics data
        HMM_MODEL.PLACE_M(:) = 0.5; % set everything to uniform probabilities
        for i=HMMParamsStruct.ChromosomesToRun
            cur_len = length(LDStruct.LD{i}.PairMat(:,1));   % Now transfer from joint to conditional probabilities
            HMM_MODEL.CHR_STARTS(i) = cur_ind; HMM_MODEL.CHR_ENDS(i) = cur_ind+cur_len;
            cur_ind = cur_ind+cur_len+1; % Take extra one, skip ..
        end

    case 1 % use single LD statsitics data
        for i=HMMParamsStruct.ChromosomesToRun
            cur_len = length(LDStruct.LD{i}.PairMat(:,1));   % Now transfer from joint to conditional probabilities
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,1) = 1-LDStruct.LD{i}.FreqVec(2:end);
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,2) = 1-HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,1);  % Switch 2 and 3 for the C function
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,3) = 1-LDStruct.LD{i}.FreqVec(2:end);
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,4) = 1-HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,3);
            HMM_MODEL.PLACE_M(cur_ind+cur_len,:) = [0.5 0.5 0.5 0.5]; % These are the 'no linkage' probs. between chromosomes
            HMM_MODEL.CHR_STARTS(i) = cur_ind; HMM_MODEL.CHR_ENDS(i) = cur_ind+cur_len;
            cur_ind = cur_ind+cur_len+1; % Take extra one, skip ..
        end

    case 2 % use 'full' LD pairwise statsitics data
        for i=HMMParamsStruct.ChromosomesToRun
            cur_len = length(LDStruct.LD{i}.PairMat(:,1));   % Now transfer from joint to conditional probabilities
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,1) = (LDStruct.LD{i}.PairMat(:,1)+HMMParamsStruct.derich) ./ (LDStruct.LD{i}.PairMat(:,1)+LDStruct.LD{i}.PairMat(:,2)+2*HMMParamsStruct.derich);
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,2) = 1-HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,1);  % Switch 2 and 3 for the C function
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,3) = (LDStruct.LD{i}.PairMat(:,3)+HMMParamsStruct.derich) ./ (LDStruct.LD{i}.PairMat(:,3)+LDStruct.LD{i}.PairMat(:,4)+2*HMMParamsStruct.derich);
            HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,4) = 1-HMM_MODEL.PLACE_M(cur_ind:cur_ind+cur_len-1,3);
            HMM_MODEL.PLACE_M(cur_ind+cur_len,:) = [0.5 0.5 0.5 0.5]; % These are the 'no linkage' probs. between chromosomes
            HMM_MODEL.CHR_STARTS(i) = cur_ind; HMM_MODEL.CHR_ENDS(i) = cur_ind+cur_len;
            cur_ind = cur_ind+cur_len+1; % Take extra one, skip ..
        end
end

