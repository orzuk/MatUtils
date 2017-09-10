% Calculate the pairwise similarity matrix between all pairs
% of pwms where one is from PWMS1 and the other from PWMS2.
% The sizes of the PWMS in PWMS2 can be different
% (but still, all the matrices in PWMS1 are assumed to be of the same
% length, and the same is true for all the matrices in PWMS2)
% One can also choose the metric to use. currently supported are Euclidian
% metric and LogLikelihood similarity,
% which gives more weight to the more informative positions in the matrix.
% The PWMS are given in the following format: Matrices of 4*L. So we need
% to convert them to 'vectors'
% The similarity is used for example for the Affinity Propagation
% clustering algorithm by Frey and Dueck
% max_shift is the maximal shift we move the pwms with respect to each
% other - we choose the shift which gives the best similarity
%
% Input:
% PWMS1 - first set of PWMS
% PWMS2 - second set of PWMS
% metric - the similarit metric to use
% rev_comp_flag - flag saying if to account also for reverse complement when doing comparison (optional)
% align_flag - flag saying if when aligning matrices you take only the overlap, or the union
% min_overlap - minimum number of nucleotide overlap between two matrices (default is one)
%
% Output:
% S - the scores of best matches to pwms
% I - the indices in the sequences giving the best scores, when computing the edit distance (should work now)
% J - orientation of 2nd pwms with respect to first pwms (when doing edit distance). 0 - means the same
%
function [S I J] = TwoPWMsToSimilarityMatrix(PWMS1, PWMS2, ...
    metric, rev_comp_flag, align_flag, min_overlap, varargin)

% NEW! this works also for cell-arrays, representing pwms of different sizes.
% The computations are, however, much slower, so its best to do this only
% for a small (< ~ 5000(?)) number of elements.
% We assume that if PWMS1 is a cell-array, then so is PWMS2.
Assign24MammalsGlobalConstants;
tolerance = 0.01; % make simialr pwms the same

if(~exist('rev_comp_flag', 'var')) % default is to try both orientations !!
    rev_comp_flag = 1;
end
if(~exist('align_flag', 'var'))
    align_flag = [];
end
if(~exist('min_overlap', 'var'))
    min_overlap = 1;
end

if(rev_comp_flag) % Need to flip only one set (not both)
    [S1 I1] = TwoPWMsToSimilarityMatrix(PWMS1, PWMS2, ...
        metric, 0, align_flag, min_overlap);
    [S2 I2] = TwoPWMsToSimilarityMatrix(PWMS1, pwmrcomplement(PWMS2), ...
        metric, 0, align_flag, min_overlap);
    [S J] =  max_with_inds(S1, S2);
    I = I1 .* (1-J) + I2 .* J; 
    return;     
%    all_pwms2 = [all_pwms2 pwmrcomplement(all_pwms2)]; %    all_pwms1 = [all_pwms1 pwmrcomplement(all_pwms1)];
end

if(iscell(PWMS1)) % New implementation: First get all columns
    [all_pwms1 all_cum_lens1] = cell2vec(PWMS1); num_pwms1 = length(PWMS1);
    all_pwms_lens1 = length_cell(PWMS1);
end
if(iscell(PWMS2))
    [all_pwms2 all_cum_lens2] = cell2vec(PWMS2); num_pwms2 = length(PWMS2);
    all_pwms_lens2 = length_cell(PWMS2); max_len2 = max(all_pwms_lens2);
    all_pwms_inds2 = cell(max_len2,1);
    for i=1:max_len2
        all_pwms_inds2{i} = find(all_pwms_lens2 >= i); 
    end
end

[unique_all_pwms1 u1 unique_inds1] = unique_with_tol(single(all_pwms1), tolerance);
[unique_all_pwms2 u2 unique_inds2] = unique_with_tol(single(all_pwms2), tolerance);

% pwms_inds2 = zeros(num_pwms2, max_len2); % Get the indicies of 2nd pwms
pwms_inds1 = vec2cell(unique_inds1, all_cum_lens1);
pwms_inds2 = vec2cell(unique_inds2, all_cum_lens2); % Get indices as cell
num_col1 = size(unique_all_pwms1,2);
num_col2 = size(unique_all_pwms2,2);

pwms_inds_shift_mat2 = zeros(num_pwms2, max_len2)+num_col2+1; % build the shift-index matrix
for i=1:num_pwms2
    pwms_inds_shift_mat2(i,1:all_pwms_lens2(i)) = pwms_inds2{i};
end
sim_mat = zeros(num_col1, num_col2, 'single'); % Get a similarity matrix for ALL columns and compute metric

switch metric % Prepare column-wise small distance matrix 
    case LOGLIKE % P_i * log(Q_i)
        sim_mat = unique_all_pwms1' * log(derich_correct(unique_all_pwms2, epsilon, 4));
    case KULLBACK_LEIBLER % P_i * log(P_i / Q_i)
        sim_mat = -unique_all_pwms1' * log(derich_correct(unique_all_pwms2, epsilon, 4)) + ...
            repmat( sum( unique_all_pwms1 .* log(derich_correct(unique_all_pwms2, epsilon, 4)) ), 1, num_col2 );
    case EUCLIDIAN % (P_i - Q_i)^2
        sim_mat = -two_vecs_to_distance(unique_all_pwms1', unique_all_pwms2').^2; 
    case DOTPROD % P_i*Q_i - 1/4 (normalization) 
        sim_mat = unique_all_pwms1' * unique_all_pwms2 - 0.25;
    case PEARSON % [\sum( P_i Q_i) - n*X*Y] / (n-1) S_x S_y  same as dot-prod but need to correct later ..
        sim_mat = unique_all_pwms1' * unique_all_pwms2;
end
sim_mat = [sim_mat zeros(num_col1,1)];

S = zeros(num_pwms1, num_pwms2); I=S; all_cum_lens1 = [0 vec2row(all_cum_lens1)];

for i=1:num_pwms1 % loop on first pwms
    pwm_inds = unique_inds1( all_cum_lens1(i)+1:all_cum_lens1(i+1) );
    pwm_mat = sim_mat(pwm_inds, :);
    
    cur_max_len = min(max_len2, all_pwms_lens1(i));
    for j=1:cur_max_len% all_pwms_lens1(i)
        S(i,all_pwms_inds2{j}) = S(i,all_pwms_inds2{j}) + pwm_mat(j, pwms_inds_shift_mat2(all_pwms_inds2{j},j) ); % Do first shift (0)
%        S(i,:) = S(i,:) + pwm_mat(j, pwms_inds_shift_mat2(:,j) ); % Do first shift (0)
    end
    I(i,:) = 0; % zero shift
    max_shift = all_pwms_lens1(i)-min_overlap; min_shift = -max_len2+min_overlap; % Temp .. 
    for shift = min_shift:max_shift % loop on both positive and negative shifts
        tmp_S = zeros(1, num_pwms2)-999; tmp_S(all_pwms_inds2{max(1,shift+1)-shift}) = 0;
        for j=max(1,shift+1):min(all_pwms_lens1(i), max_len2+shift)
%            tmp_S = tmp_S + pwm_mat(j, pwms_inds_shift_mat2(:,j-shift) ); % Do other shift (0)
            tmp_S(all_pwms_inds2{j-shift}) = tmp_S(all_pwms_inds2{j-shift}) + ...
                pwm_mat(j, pwms_inds_shift_mat2(all_pwms_inds2{j-shift},j-shift) ); % Do other shift (0)
        end
        [S(i,:) S_inds] = max_with_inds(S(i,:), tmp_S);
        I(i, find(S_inds)) = shift;
    end
end

switch metric % Final metric-specific corrections
    case EUCLIDIAN
        S = -sqrt(-S); 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here: Old Stuff ...
do_old = 0;
if(do_old)
    if(iscell(PWMS1))
        S = zeros(length(PWMS1), length(PWMS2)); I=S;
        len_vec = zeros(length(PWMS1), 1);
        for i=1:length(PWMS1)
            len_vec(i) = length(PWMS1{i});
        end
        unique_len_vec = unique(len_vec);
        
        for i=1:length(unique_len_vec)
            sprintf('doing len # %d out of %d', i, length(unique_len_vec))
            cur_len_inds = find(len_vec == unique_len_vec(i));
            cur_pwms1 = zeros(4, unique_len_vec(i), length(cur_len_inds));
            for j=1:length(cur_len_inds)
                cur_pwms1(:,:,j) = PWMS1{cur_len_inds(j)};
            end
            [S(cur_len_inds,:) I(cur_len_inds,:)] = ...
                TwoPWMsToSimilarityMatrix(pwms_to_vecs(cur_pwms1), PWMS2, metric, rev_comp_flag, align_flag, min_overlap);
        end
    else    % PWM1 is not cell
        if(rev_comp_flag)
            PWMS2_RC = PWMS2; n2 = length(PWMS2);
            for i=1:n2
                PWMS2_RC{i+n2} = pwmrcomplement(PWMS2{i});
            end
            [S I] = TwoPWMsToSimilarityMatrix(PWMS1, PWMS2_RC, metric, 0, align_flag, min_overlap);
            inds_mat = sign(S(:,1:n2) - S(:,n2+1:2*n2));
            S = max(S(:,1:n2), S(:,n2+1:2*n2));
            I = I(:,1:n2) .* inds_mat + I(:,n2+1:2*n2) .* (1-inds_mat);
        else
            if(length(size(PWMS1)) == 3)
                V1 = pwms_to_vecs(PWMS1); % V2 = pwms_to_vecs(PWMS2);
                L = size(PWMS1, 2);
            else
                V1 = PWMS1; L = size(PWMS1, 1) / 4;
            end
            
            num_pwms1 = size(V1, 2);
            num_pwms2 = length(PWMS2);
            S = zeros(num_pwms1, num_pwms2); I=S;
            switch metric
                % %             case EUCLIDIAN % euclidian - NOT WORKING YET !!!!
                % %                 % Finally, computing the distance function for the matrix is really simple
                % %                 num_pwms1 = size(V1, 2);
                % %                 num_pwms2 = length(PWMS2);
                % %                 S = -vecs_to_distance(V1');
                case LOGLIKE  % Log-Likelihhod \sum_i P_i log Q_i (not symmetric)
                    epsilon = 0.00000000001; V1 = derich_correct(V1, epsilon, 4); %%   (V1 + epsilon) ./ (1+4*epsilon);
                    for i=1:num_pwms2 % loop on 2nd pwms
                        L2 = size(PWMS2{i}, 2);
                        V2 = pwms_to_vecs(PWMS2{i});       %%            V2 = reshape(PWMS2{i}, L2*4, 1);
                        V2 = derich_correct(V2, epsilon, 4); %%                    V2 = (V2 + epsilon) ./ (1+4*epsilon);
                        if(L2 <= L)
                            S(:,i) = V1(1:4*L2,:)' * log(V2);
                            max_shift = L-L2;
                            for s=1:max_shift
                                I( find(S(:,i) < V1(s*4+1:s*4+4*L2,:)' * log(V2)), i ) = s;
                                S(:,i) = max( S(:,i), V1(s*4+1:s*4+4*L2,:)' * log(V2) );
                            end
                        else % here L2 < L
                            S(:,i) = V1' * log(V2(1:4*L,:));
                            max_shift = L-L2;
                            for s=1:max_shift
                                I( find(S(:,i) < V1(s*4+1:s*4+4*L2,:)' * log(V2)), i ) = s;
                                S(:,i) = max( S(:,i), V1' * log(V2(s*4+1:s*4+4*L,:)) );
                            end
                        end
                    end
                case KULLBACK_LEIBLER
                    
                case {EUCLIDIAN,PEARSON,DOTPROD} % 0 - (minus) euclidian distance 2 - pearson correlation, 3 - simple dot-product
                    if(metric == 2)
                        V1 = V1 - repmat(mean(V1), size(V1, 1), 1);
                        V1 = V1 ./ repmat(sqrt(sum(V1.^2)), size(V1,1), 1);
                    end
                    metric_normalise = max(metric,2)-2; % if we take dot-product we need to normalize by length
                    for i=1:num_pwms2
                        %i_is = i;
                        L2 = size(PWMS2{i}, 2);
                        V2 = reshape(PWMS2{i}, L2*4, 1); % this is just like pwms_to_vecs
                        if(metric == 2)
                            V2 = V2 - mean(V2);
                            V2 = V2 ./ sqrt(sum(V2.^2)); % normalized
                        end
                        %max_shift = abs(L2-L); % allow shift
                        
                        [W1 W2] = TwoPWMsToShiftComplete(V1, V2, 0, align_flag); % either intersection or union
                        if(L2 <= L) % initilize
                            if(metric == EUCLIDIAN) % Use normalized euclidian distance
                                S(:,i) = -sum( (W1(1:4*L2,:) - W2).^2 ) / L2;
                            else
                                S(:,i) = W1(1:4*L2,:)' * W2 - metric_normalise*L2/4;
                            end
                        else
                            if(metric == EUCLIDIAN)
                                S(:,i) = -sum( (V1 - V2(1:4*L,:)).^2 ) / L;
                            else
                                S(:,i) = V1' * V2(1:4*L,:) - metric_normalise*L/4;
                            end
                        end
                        %min_L = min(L,L2); % take the minimal length
                        for s=1:L-min_overlap % slide the second
                            [W1 W2] = TwoPWMsToShiftComplete(V1, V2, s, align_flag); % either intersection or union
                            
                            if(metric == EUCLIDIAN)
                                cur_shift_score = -sum( (W1(s*4+1:4*min(L,L2+s),:) - W2(1:4*min(L2,L-s),:)).^2 ) / min(L2,L-s);
                            else
                                cur_shift_score = W1(s*4+1:4*min(L,L2+s),:)' * W2(1:4*min(L2,L-s),:) - ...
                                    metric_normalise*min(L2,L-s)/4
                            end
                            I( find(S(:,i) < cur_shift_score), i ) = s;
                            S(:,i) = max( S(:,i), cur_shift_score );
                        end
                        for s=1:L2-min_overlap % slide the first
                            [W1 W2] = TwoPWMsToShiftComplete(V1, V2, -s, align_flag); % either intersection or union
                            if(metric == EUCLIDIAN)
                                cur_shift_score = -sum( (W1(1:4*min(L,L2-s),:) - W2(s*4+1:4*min(L2,L+s),:)).^2 ) / min(L2,L+s);
                            else
                                cur_shift_score = W1(1:4*min(L,L2-s),:)' * W2(s*4+1:4*min(L2,L+s),:) - ...
                                    metric_normalise*min(L,L2-s)/4
                            end
                            I( find(S(:,i) < cur_shift_score), i ) = -s;
                            S(:,i) = max( S(:,i), cur_shift_score );
                        end
                    end
                    if(metric == EUCLIDIAN) % take the sqrt
                        S = -sqrt(abs(S));
                    end
            end % switch metric
        end % if rev comp flag
    end
end % old stuff

% This subfunction generates the appropriate vectors, given a shift
% by completing stuff not aligned to the background distribution (0.25 each)
function [W1 W2] = TwoPWMsToShiftComplete(V1, V2, shift, align_flag)
[m1 n1] = size(V1); [m2 n2] = size(V2);

end_ind = max(L, L2+shift); % total aligned set
start_ind = min(1, shift+1);
total_L = end_ind - start_ind+1;

if(align_flag) % take union
    W1 = [0.25.*ones(m1, 4*(1-start_ind))  V1 0.25.*ones(m1, 4*(end_ind-L))];
    W2 = [0.25.*ones(m2, 4*(1-start_ind))  V2 0.25.*ones(m2, 4*(end_ind-L2))];
else % take intersection - NEED TO FIX THIS!!!
    W1 = V1;
    W2 = V2;
end

% switch metric
%     case 1 % Log-Likelihhod \sum_i P_i log Q_i (not symmetric)
%
%
%
%     case 2
%
%
%     case 3
%
% end
%S=[];I=[];
