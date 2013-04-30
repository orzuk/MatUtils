% Compute all inference functions obtained by a finite-length binary HMM  
% given an observed set. 
% Compute all the possible maximal sequences
% for a whole range of parameters eps_vec, p_vec.
% What do we maximize? according to the flag :
% If flag = find_Y   (0) : find Y maximizing P(Y).
% If flag = finx_X   (1)  : find X maximizing P(X,Y) for a given Y
% If flag = find_XY (2):  find (X,Y) maximizing P(X,Y)
%
% inf_flag - which inference function to use: 0 - MLE 1 - Marginal 
function [min_seqs, min_vals, max_seqs, max_vals, NF] = HMPInferenceFunction(eps_vec, p_vec, N, seqs_flag, Y, inf_flag, varargin)

NF = [];

min_seqs = zeros(length(p_vec), length(eps_vec)); 
min_vals = zeros(length(p_vec), length(eps_vec)); 
max_seqs = zeros(length(p_vec), length(eps_vec)); 
max_vals = zeros(length(p_vec), length(eps_vec)); 
TOL = 0.0000000001;

MLE = 0; MARGINAL = 1;
if(~exist('inf_flag', 'var'))
    inf_flag = MLE; 
end


% we loop over p and send to the P_Y function
switch seqs_flag 
    case 0 % maximize over all Y's
        for p=1:length(p_vec)
            p_ind_is = p
            cur_p_vec = p_vec(p) + zeros(size(eps_vec, 1), size(eps_vec, 2));
            P_Y = HMP_ProbY(eps_vec, cur_p_vec, N);
            [max_vals(p,:) max_seqs(p,:)] = max(P_Y(:,1:2^(N-1)),[],2);     
            [min_vals(p,:) min_seqs(p,:)] = min(P_Y(:,1:2^(N-1)),[],2);     
        end
    case 1 % maximize over all X's for a given Y
        for p=1:length(p_vec)
            p_is = p
            cur_p_vec = p_vec(p) + zeros(size(eps_vec, 1), size(eps_vec, 2));
            if(inf_flag == MLE)
                P_X_given_Y = HMP_ProbX_given_Y(eps_vec, cur_p_vec, N, Y);
                [max_vals(p,:) max_seqs(p,:)] = max_with_tol(P_X_given_Y(:,1:2^(N-1)),2, TOL);   % let tolerance when maximising
                [min_vals(p,:) min_seqs(p,:)] = min(P_X_given_Y(:,1:2^(N-1)),[],2);
            else % here marginal
                P_X_given_Y_marginal = HMP_ProbX_given_Y_marginal(eps_vec, cur_p_vec, N, Y);
                for j=1:1 % N
                    [tmp_max_vals tmp_max_seqs] = max_with_tol(P_X_given_Y_marginal(:,j*2-1:j*2),2, TOL);   % let tolerance when maximising
                    max_vals(p,:) = max_vals(p,:) + vec2row(tmp_max_vals);
                    max_seqs(p,:) = max_seqs(p,:) + 2^(j-1) * vec2row(tmp_max_seqs);
                end
            end
        end
    case 2 % maximize over all pairs (X,Y)
        
        
    case 3 % Here be a bit more series and start looking at convex hulls etc.
        NP = NewtownPolytope(N,Y); % compute the Newton Polytope
        NF = NormalFan(NP); % Compute the normal fan of the Newton Polytope
        
        plot_polytope=1;
        if(plot_polytope)
            subplot(2^floor(N/2),2^ceil(N/2),Y+1);  % figure;
            plot(NP(:,1), NP(:,2), 'linewidth', 3); % Plot the Newton Polytope
            axis( [min(NP(:,1))-1 max(NP(:,1))+1 min(NP(:,2))-1 max(NP(:,2))+1]);
            title(['Newton Polytope of Y = ' num2str(Y)]);
        else % here plot fan
            subplot(2^floor(N/2),2^ceil(N/2),Y+1);  hold on; % figure;
            for i=1:size(NF,1)
                line(  [0 NF(i,1)], [0,NF(i,2)], 'linewidth', 3); % Plot the Newton Polytope
            end
            circle([0,0], 1, 500, '.');
%            axis( [min(NP(:,1))-1 max(NP(:,1))+1 min(NP(:,2))-1 max(NP(:,2))+1]);
            title(['Normal Fan of Y = ' num2str(Y)]);
        end
            
        
        
end

% boundary_inds = find(diff(min_seqs));

function [max_val max_ind] = max_with_tol(V, dim, TOL, varargin)

if(~exist('dim', 'var'))
    dim = 1;
end
if(~exist('TOL', 'var'))
    TOL = 0.0000000001;
end

[max_val max_ind] = max(V,[],dim); % perform regular max

for i=1:size(V,1)
    max_ind(i) = find(V(i,:) > max_val(i) - TOL , 1); % take the first index, lexicographicaly % repmat(max_val,1,size(V,2)) - TOL
    max_val(i) = V(i,max_ind(i));
end
%max_val = V(max_ind);



