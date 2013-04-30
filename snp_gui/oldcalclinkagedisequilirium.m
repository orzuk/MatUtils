% Calculate the Linkage Disequilibrium values of all pairs
% FreqData - The (singleton) SNPs frequencies - currently not used but
% calculated directly from the data
% DataMat - The data matrix
% OffspringsInds - The childs we ignore in the statistics (can be empty)
%
function [LD_mat LocalFreqVec PairwiseFreqs ] = ...
    CalcLinkageDisequilirium(FreqVec, DataMat, BadCallsMat, OffspringsInds, LD_flag, use_hetro_flag, Calc_LD_mat_flag)
TOL = 0.00000000000000000001;

LD_mat = {};

N=length(FreqVec); % Number of SNPs

% Try a faster way ...
if(Calc_LD_mat_flag)
    LD_mat = zeros(N);
    LD_BadCounts = zeros(N);
    LocalFreq_mat1 = zeros(N); LocalFreq_mat2 = zeros(N);% We need a different FREQUENCY for each pair!!!!
end

% Get all the pairwise frequencies

PairwiseFreqs = {};
PairwiseFreqs{1} = zeros(N-1,1); % 00
PairwiseFreqs{2} = zeros(N-1,1); % 01
PairwiseFreqs{3} = zeros(N-1,1); % 10
PairwiseFreqs{4} = zeros(N-1,1); % 11
PairwiseBadCounts = zeros(N-1,1); % How many are baddd


% First calculate the frequency vector again!!!
LocalFreqVec = zeros(N,1);
BadCounts = zeros(N,1);

if(use_hetro_flag == 1) % Use all data assuming it is paired correctly
    for j=1:6  % Each bit represent a person and his allele
        for b=1:30
            if( isempty(find(OffspringsInds == mod(j-1,3)*30+ b)) )
                A = 1-bitget(DataMat(:,j), b); % Major allele: A
                BadCallsVec = bitget(BadCallsMat(:,j), b);
                LocalFreqVec = LocalFreqVec + A .* (1-BadCallsVec);
                BadCounts = BadCounts + BadCallsVec;
                if(Calc_LD_mat_flag)
                    LocalFreq_mat1 = LocalFreq_mat1 + ((1-BadCallsVec)*(1-BadCallsVec)') .* repmat(A,1,N);
                    LocalFreq_mat2 = LocalFreq_mat2 + ((1-BadCallsVec)*(1-BadCallsVec)') .* repmat(A',N,1);
                    LD_BadCounts = LD_BadCounts + 1 - (1-BadCallsVec)*(1-BadCallsVec)'; % Add one to the counts
                    LD_mat = LD_mat + (A .* (1-BadCallsVec)) * (A .* (1-BadCallsVec))';
                end
                PairwiseBadCounts = PairwiseBadCounts + 1 - (1-BadCallsVec(1:end-1)).*(1-BadCallsVec(2:end));
                tmp_A = A .* (1-BadCallsVec); tmp_one_minus_A = (1-A) .* (1-BadCallsVec);
                PairwiseFreqs{1} = PairwiseFreqs{1} + tmp_one_minus_A(1:end-1) .* tmp_one_minus_A(2:end);
                PairwiseFreqs{2} = PairwiseFreqs{2} + tmp_one_minus_A(1:end-1) .* tmp_A(2:end);
                PairwiseFreqs{3} = PairwiseFreqs{3} + tmp_A(1:end-1) .* tmp_one_minus_A(2:end);
                PairwiseFreqs{4} = PairwiseFreqs{4} + tmp_A(1:end-1) .* tmp_A(2:end);
            end
        end
    end

    LocalFreqVec = LocalFreqVec./(120.0-BadCounts); % 180.0
    if(Calc_LD_mat_flag)
        LocalFreq_mat1 = LocalFreq_mat1 ./ (120.0-LD_BadCounts);
        LocalFreq_mat2 = LocalFreq_mat2 ./ (120.0-LD_BadCounts);
        LD_mat = LocalFreq_mat1 .* LocalFreq_mat2 - LD_mat./(120.0-LD_BadCounts); %180.0
    end

    %    dddd = diag(LD_BadCounts(2:end,1:end-1))
    PairwiseFreqs{1} = PairwiseFreqs{1} ./ (120.0-PairwiseBadCounts);
    PairwiseFreqs{2} = PairwiseFreqs{2} ./ (120.0-PairwiseBadCounts);
    PairwiseFreqs{3} = PairwiseFreqs{3} ./ (120.0-PairwiseBadCounts);
    PairwiseFreqs{4} = PairwiseFreqs{4} ./ (120.0-PairwiseBadCounts);

    % Normalize by st.d.
    if(Calc_LD_mat_flag)
        LD_mat = LD_mat ./ sqrt(max(LocalFreq_mat1,TOL));
        LD_mat = LD_mat ./ sqrt(max(1-LocalFreq_mat1,TOL));
        LD_mat = LD_mat ./ sqrt(max(LocalFreq_mat2,TOL));
        LD_mat = LD_mat ./ sqrt(max(1-LocalFreq_mat2,TOL));
        mean_BAD_counts = mean(mean(LD_BadCounts))
    end

    % % %     % Normalize by st.d.
    % % %     LD_mat = LD_mat ./ repmat(sqrt(max(LocalFreqVec,TOL)),1,N);
    % % %     LD_mat = LD_mat ./ repmat(sqrt(max(1-LocalFreqVec,TOL)),1,N);
    % % %     LD_mat = LD_mat ./ repmat(sqrt(max(LocalFreqVec,TOL)'),N,1);
    % % %     LD_mat = LD_mat ./ repmat(sqrt(max(1-LocalFreqVec,TOL)'),N,1);

else  % Throw away data in which both SNPs are heterozygoces
    Freq_hetro_counts = zeros(N,1);
    LD_hetro_counts = zeros(N); % A matrix counting how many times we've got an hetrozygoce couple
    for j=1:3  % We take a person's two allels together
        for b=1:30
            if( isempty(find(OffspringsInds == mod(j-1,3)*30+ b)) )
                A = 1-bitget(DataMat(:,j), b); % first allele
                B = 1-bitget(DataMat(:,j+3), b); % second allele
                BadCallsVecA = bitget(BadCallsMat(:,j), b);
                BadCallsVecB = bitget(BadCallsMat(:,j+3), b);
                maxi = max((BadCallsVecA-BadCallsVecB).^2)
                BadCallsVec = 1 .* bitor(BadCallsVecA, BadCallsVecB); % Throw away even if one N error!
                LocalFreqVec = LocalFreqVec + A+B; % This is not affected
                hetro_vec = zeros(N,1); hetro_vec(find(A ~=B))=1;  % Find the hetrozygotes SNPs for this person
                LocalFreq_mat1 = LocalFreq_mat1 + (1-hetro_vec*hetro_vec').* ...
                    ((1-BadCallsVec)*(1-BadCallsVec)') .* repmat(A+B,1,N);
                LocalFreq_mat2 = LocalFreq_mat2 + (1-hetro_vec*hetro_vec').* ...
                    ((1-BadCallsVec)*(1-BadCallsVec)') .* repmat(A'+B',N,1);
                %%%%                LocalFreq_mat1 = LocalFreq_mat1 + ((1-BadCallsVec)*(1-BadCallsVec)') .* repmat(A+B,1,N);
                %%%%                LocalFreq_mat2 = LocalFreq_mat2 + ((1-BadCallsVec)*(1-BadCallsVec)') .* repmat(A'+B',N,1);
                if(Calc_LD_mat_flag)
                    LD_mat = LD_mat + (1-hetro_vec*hetro_vec') .* ((1-BadCallsVec)*(1-BadCallsVec)') .* (A*A'+B*B');
                end
                %%%%                    0.5 .* (hetro_vec*hetro_vec') .* ((1-BadCallsVec)*(1-BadCallsVec)'); % The hetros gave count of 2 which we subtract back
                LD_hetro_counts = LD_hetro_counts + 2 .* bitor( (hetro_vec*hetro_vec'), ...    % Add one to the counts
                    (1 - (1-BadCallsVec)*(1-BadCallsVec)') ); % Add one to the counts
            end
        end
    end
    LocalFreqVec = LocalFreqVec./120.0; % 180.0 % Not affected
    if(Calc_LD_mat_flag)
        LocalFreq_mat1 = LocalFreq_mat1 ./ (120.0-LD_hetro_counts);
        LocalFreq_mat2 = LocalFreq_mat2 ./ (120.0-LD_hetro_counts);
        LD_mat = LocalFreq_mat1 .* LocalFreq_mat2 - LD_mat./(120.0-LD_hetro_counts);
    end

    % Normalize by st.d.
    if(Calc_LD_mat_flag)
        LD_mat = LD_mat ./ sqrt(max(LocalFreq_mat1,TOL));
        LD_mat = LD_mat ./ sqrt(max(1-LocalFreq_mat1,TOL));
        LD_mat = LD_mat ./ sqrt(max(LocalFreq_mat2,TOL));
        LD_mat = LD_mat ./ sqrt(max(1-LocalFreq_mat2,TOL));
    end

    %    LD_hetro_counts_is = LD_hetro_counts
    LD_BadCounts = LD_hetro_counts; %%%% OUTPUT THE BAD COUNTS (HETRO)
    mean_hetro_counts = mean(mean(LD_hetro_counts))
end

if(Calc_LD_mat_flag)
    for i=1:N
        LD_mat(i,1:i-1)=0;
    end
end
