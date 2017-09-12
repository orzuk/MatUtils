% Draw at random the vertex which is chosen. Do it according to the poisson
% rate
function [flipped_ind flip_time lambda_total] =randomize_vertex_uniform(X, lambda_vec)

N = length(X); % # of vertices
N_on = sum(X); % Initial # of infected vertices

X_off = find(X==0);

 lambda_total= sum(lambda_vec(X_off));
length(X_off)
 
r = rand(1)*lambda_total; % randomize

flipped_ind = ceil(r); 

% search which vertex we should flip. Do a binary search
% % % if(r < lambda_off_cumsum_vec(1))
% % %     flipped_ind=1;
% % % else
% % %     left_ind=1; right_ind=length(lambda_off_cumsum_vec);
% % % 
% % %     while(right_ind-left_ind>1)
% % %         mid_ind = round((left_ind+right_ind)/2);
% % %         if(r < lambda_off_cumsum_vec(mid_ind))
% % %             right_ind=mid_ind;
% % %         else
% % %             left_ind=mid_ind;
% % %         end
% % %      
% % %     end
% % %        
% % %     flipped_ind = right_ind; 
% % % end


% Transform back to the original indexes
flipped_ind = X_off(flipped_ind);

% Randomize also the time until flip
flip_time = exprnd(1/lambda_total);