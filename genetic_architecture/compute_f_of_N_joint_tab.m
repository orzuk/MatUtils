% Compute joint probability distribution of one varialbe and the entire out
% of a binomial 
function joint_tab = compute_f_of_N_joint_tab(N, K, p)
    joint_tab = zeros(2,2); 

    mu = 1-binocdf(K-1,N,p); % probability at least K are 'on'
    joint_tab(2,2) = p* (1-binocdf(K-2,N-1,p)); % both vars are set to one
    joint_tab(1,2) = p-joint_tab(2,2);
    joint_tab(2,1) = mu-joint_tab(2,2);
    joint_tab(1,1) = 1-joint_tab(1,2)-joint_tab(2,1)-joint_tab(2,2);
end
