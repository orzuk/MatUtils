% Calculate the pairwise similarity matrix between all pairs
% of pwms where one is from PWMS 1 and the other from PWMS2.
% The difference here is that the sizes of the PWMS in PWMS2 can be different
% One can choose the metric to use. currently supported are Euclidian
% metric and LogLikelihood similarity,
% which gives more weight to the more informative positions in the matrix.
% The PWMS are given in the following format: Matrices of 4*L. So we need
% to convert them to 'vectors'
% The similarity is used for example for the Affinity Propagation
% clustering algorithm by Frey and Dueck
% max_shift is the maximal shift we move the pwms with respect to each
% other - we choose the shift which gives the best similarity 
% We give all shifts score - we assume we just do here one sequence 
% (so one PWMS1) but the number of pwms (PWMS2) can be anything
function [S I] = TwoPWMsToSimilarityVector(PWMS1, PWMS2, metric)

if(length(size(PWMS1)) == 3)
    V1 = pwms_to_vecs(PWMS1); % V2 = pwms_to_vecs(PWMS2);
    L = size(PWMS1, 2); 
else
    V1 = PWMS1; L = size(PWMS1, 1) / 4; 
end

if(metric == 0) % euclidian - NOT WORKING YET !!!!
    % Finally, computing the distance function for the matrix is really
    % simple
    S = -vecs_to_distance(V');
else  % Log-Likelihhod
    num_pwms1 = size(V1, 2)
    num_pwms2 = length(PWMS2)
    S = zeros(num_pwms1, num_pwms2); I=S;
    for i=1:num_pwms2
        i_is = i;
        L2 = size(PWMS2{i}, 2);
        V2 = reshape(PWMS2{i}, L2*4, 1);
        %%        V2 = pwms_to_vecs(PWMS2{i});
        epsilon = 0.00000000001;
        V1 = (V1 + epsilon) ./ (1+4*epsilon); V2 = (V2 + epsilon) ./ (1+4*epsilon);
        if(L2 <= L)
            max_shift = L-L2; S = zeros(1, max_shift); I = 1:max_shift+1;
            S(1) = V1(1:4*L2,:)' * log(V2);
            for s=1:max_shift
                S(s+1) = max( V1(s*4+1:s*4+4*L2,:)' * log(V2) );
            end
        else % here L2 < L
            max_shift = L-L2; S = zeros(1, max_shift); I = 1:max_shift+1;
            S(1) = V1' * log(V2(1:4*L,:));
            for s=1:max_shift
                S(s+1) = max( V1' * log(V2(s*4+1:s*4+4*L,:)) );
            end
        end
    end
end





