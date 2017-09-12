function scores = my_BIC_score_all_dags(data_counts, num_samples, dags)
% my_BIC_SCORE_ALL_DAGS Compute the mdl BIC score of many DAGs
% given counts from data, so we can save computation by letting each 
% part of the score be computed seperately.
% score = score_dags(data, ns, dags, varargin)
%
% data_count = value of all nodes combinations
% num_samples = number of samples used in this sample
% dags{g} is the g'th dag
% score(g) is the score of the i'th dag
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%

TOL = 0.00000000000000000000000001;

NG = length(dags); % number of graphs
num_iters = size(data_counts,1); % number of different sampling vectors 
n = length(dags{1}); % number of nodes


scores = zeros(NG, num_iters);

% The distribution of the counts. Assume at least they all have the same
% number of samples 
P = max((data_counts ./ num_samples)', TOL);



% Loop over the possible sets of nodes, 
% and compute the entropy of each ..
size(P,2)
H_P = zeros(2^n, size(P,2));

ham_vec = hamming_weight(0:2^n-1);

for set_i=0:2^n-1    
    set_comp_i = 2^n-1-set_i;
%     ind_i = find(bitget(set_i, 1:n));
%     ind_comp_i = setdiff(1:n, ind_i);
%     
%     bbbbbbbb = bitxor(set_i, set_comp_i)
    SUM_P_x = zeros(1, size(P,2));
    
    % Loop over the X's
    for j=0:2^ham_vec(set_i+1)-1
        P_x = zeros(1, size(P,2));
        
        % Loop over the Y's
        for k=0:2^(n-ham_vec(set_i+1))-1
%             sss_j= stretch_to_indexes(j,set_i)
%             sss_k = stretch_to_indexes(k, set_comp_i)
%             is_zero_bitand = bitand(sss_j, sss_k)
            
            cur_xy = stretch_to_indexes(j,set_i) + strech_to_indexes(k, set_comp_i) + 1;            
            
            P_x = P_x + P(cur_xy,:);
        end                
        SUM_P_x = SUM_P_x+P_x;
              
        % update the entropy 
        H_P(set_i+1,:) = H_P(set_i+1,:) + (P_x .* log(P_x));
    end            
end


% P_IS = P
% SUM_P_x
% H_P
% HP_MEAN = mean(H_P)
% HP_MEAN2 = mean(H_P, 2)

% loop over dags
for g=1:NG

    LL_G = 0;
    % loop over nodes
    for i=1:n       
        % identify the parents set
        pa_i = find(dags{g}(:,i));        
        pa_i_bin = sum(2.^(pa_i-1));
        pa_i_and_i_bin = sum(2.^ ([pa_i' i]-1)); %%%  pa_i_bin + 2^(i-1);
        
        LL_G = LL_G + H_P(pa_i_and_i_bin+1,:) - H_P(pa_i_bin+1,:);      
    end
    
    % Find the distribution closest to P on G
%    P_G = project_on_bnet(P, dags{g});
    
    % Give the likelihood score
%    LL_G = num_samples * sum(P .* log(P_G));
    
%     LL_COOL_G = num_samples .* LL_G 

    % Reduce BIC penalty term
    scores(g,:) = num_samples .* LL_G  - dim_dag(dags{g}) * log(num_samples)/2; 
    
end


