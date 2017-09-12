% Compute for many (discrete noded) bnets, the mutual information
% for every pair of variables. Currently works only for binary nodes
function bnets_MI_mats = bnets_to_pairwise_MI(bnets, way)
TOL = 0.000000000000001;

if(way == 0) % old way
    % Choose inference engine. Note : this enumerative engine is very slow, and
    % should be replaced by a different one. I tried to use jtree or var_elim,
    % but for jtree there are sets of nodes which cannot be margenalized, and
    % for var_elim, which seemed to me the best engine to use, i got an error
    % (something about no T2 variable) - could be a bug!
    nodes = length(bnets{1}.dag);
    % j=1; % Weird matlab bug


    evidence = cell(1, nodes);
    epsilon = 0.000000001;
    bnets_MI_mats = zeros(nodes, nodes, length(bnets));

    % Open all the nodes ..
    for i=1:length(bnets)
        engine = var_elim_inf_engine(bnets{i});
        engine = enumerative_inf_engine(bnets{i});

        %     engine = pearl_inf_engine(bnets{i},  'protocol', 'tree'); %%%   engine = pearl_inf_engine(bnets{i});
        %     engine = pearl_inf_engine(bnets{i},  'protocol', 'parallel'); %%%   engine = pearl_inf_engine(bnets{i});
        engine = enter_evidence(engine, evidence);

        % Get the marginal over every node
        for j=1:nodes
            m = marginal_nodes(engine,j);
            marginal_entropy(j) = entropy(m.T);
            mmm(j) = m.T(1);
            
        end
        for j=1:nodes
            for k=j+1:nodes
                mm = marginal_nodes(engine,[j,k]);
                bnets_MI_mats(j,k,i) = entropy(mm.T(:));
            end
        end
        
        bnets_MI_mats(:,:,i) = -bnets_MI_mats(:,:,i) - bnets_MI_mats(:,:,i)' + ...
            repmat(marginal_entropy, nodes,1) + repmat(marginal_entropy', 1,  nodes) - diag(marginal_entropy);
        
    end






else % my way
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Try to do the function my way: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = length(bnets{1}.dag);
    
    
    
    for i=1:length(bnets)
            % Order everything according o cannonical order
        bnets{i}.rev_order = zeros(1,n);
        for j=1:n
            bnets{i}.rev_order(bnets{i}.order(j)) = j;
        end
        
        bnets{i}.dag = bnets{i}.dag(bnets{i}.order, bnets{i}.order);
        bnets{i}.my_CPT = bnets{i}.my_CPT(bnets{i}.order,:);
    %    bnets{i}.rev_order = bnets{i}.order;
        bnets{i}.order = 1:n;

        
        % Get the parents of each variable
        T_I = find(bnets{i}.dag); T_J = ceil(T_I/n); T_I = mod(T_I-1,n)+1;
        T_Parents = zeros(1,n);
        for j=1:n
            tmp = find(bnets{i}.dag(:,j));
            if(~isempty(tmp))
                T_Parents(j) = tmp;
            end
        end

        % First get all the marginals
        marginal_probs_vec = zeros(1,n);
        marginal_probs_vec(bnets{i}.order(1)) = bnets{i}.my_CPT(bnets{i}.order(1),1);
        % Loop over all the others
        for j=2:n
            marginal_probs_vec(bnets{i}.order(j)) = marginal_probs_vec(T_Parents(bnets{i}.order(j))) * ...
                bnets{i}.my_CPT(bnets{i}.order(j),1) + ...
                (1-marginal_probs_vec(T_Parents(bnets{i}.order(j)))) * bnets{i}.my_CPT(bnets{i}.order(j),2);
        end

        joint_prob_mat = zeros(n,n,4);



        % First go over all adjacent vertices
        for j=1:n % First vertex
            for k=j+1:n % Second vertex
                prev_k = T_Parents(bnets{i}.order(k));
                if(prev_k == bnets{i}.order(j)) % Adjacent vertices
                    first_ver = min(bnets{i}.order(j), bnets{i}.order(k));
                    second_ver = max(bnets{i}.order(j), bnets{i}.order(k));

                    joint_prob_mat(first_ver,second_ver,1) = ...
                        marginal_probs_vec(bnets{i}.order(j)) * bnets{i}.my_CPT(bnets{i}.order(k),1);
                    joint_prob_mat(first_ver,second_ver,2) = ...
                        marginal_probs_vec(bnets{i}.order(j)) * (1-bnets{i}.my_CPT(bnets{i}.order(k),1));
                    joint_prob_mat(first_ver,second_ver,3) = ...
                        (1-marginal_probs_vec(bnets{i}.order(j))) * bnets{i}.my_CPT(bnets{i}.order(k),2);
                    joint_prob_mat(first_ver,second_ver,4) = ...
                        (1-marginal_probs_vec(bnets{i}.order(j))) * (1-bnets{i}.my_CPT(bnets{i}.order(k),2));
                end
            end
        end

        joint_prob_mat(:,:,1) = joint_prob_mat(:,:,1) + joint_prob_mat(:,:,1)';
        joint_prob_mat(:,:,2) = joint_prob_mat(:,:,2) + joint_prob_mat(:,:,2)';
        joint_prob_mat(:,:,3) = joint_prob_mat(:,:,3) + joint_prob_mat(:,:,3)';
        joint_prob_mat(:,:,4) = joint_prob_mat(:,:,4) + joint_prob_mat(:,:,4)';

        % Now go over non-adjacent vertices
        for j=1:n % First vertex
            for k=j+1:n % Second vertex
                prev_k = T_Parents(bnets{i}.order(k));
                if(prev_k ~= bnets{i}.order(j)) % Non-Adjacent vertices
                    % We need to decide if prev_k is before or after j  0->0->0 or 0<-0->0
                    prev_j = bnets{i}.order(j); prev_k2 = prev_k;

                    if(prev_j > prev_k) % prev_k) % SWAP if neccessary
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),1) = ...
                            joint_prob_mat(prev_j,prev_k2,1) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),1) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,3) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),3) / (1-marginal_probs_vec(prev_k));
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),2) = ...
                            joint_prob_mat(prev_j,prev_k2,1) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),2) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,3) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),4) / (1-marginal_probs_vec(prev_k));
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),3) = ...
                            joint_prob_mat(prev_j,prev_k2,2) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),1) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,4) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),3) / (1-marginal_probs_vec(prev_k));
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),4) = ...
                            joint_prob_mat(prev_j,prev_k2,2) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),2) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,4) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),4) / (1-marginal_probs_vec(prev_k));
                    else % 'Regular' ordering
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),1) = ...
                            joint_prob_mat(prev_j,prev_k2,1) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),1) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,2) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),3) / (1-marginal_probs_vec(prev_k));
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),2) = ...
                            joint_prob_mat(prev_j,prev_k2,1) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),2) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,2) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),4) / (1-marginal_probs_vec(prev_k));
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),3) = ...
                            joint_prob_mat(prev_j,prev_k2,3) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),1) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,4) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),3) / (1-marginal_probs_vec(prev_k));
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),4) = ...
                            joint_prob_mat(prev_j,prev_k2,3) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),2) / marginal_probs_vec(prev_k) + ...
                            joint_prob_mat(prev_j,prev_k2,4) * ...
                            joint_prob_mat(prev_k,bnets{i}.order(k),4) / (1-marginal_probs_vec(prev_k));
                    end

                    %                         % Swap if needed
                    if(bnets{i}.order(j) > bnets{i}.order(k))
                        tmp = joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),3);
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),3) = joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),2);
                        joint_prob_mat(bnets{i}.order(j),bnets{i}.order(k),2) = tmp;
                    end


                    joint_prob_mat(:,:,1) = max(joint_prob_mat(:,:,1),joint_prob_mat(:,:,1)');
                    joint_prob_mat(:,:,2) = max(joint_prob_mat(:,:,2),joint_prob_mat(:,:,2)');
                    joint_prob_mat(:,:,3) = max(joint_prob_mat(:,:,3),joint_prob_mat(:,:,3)');
                    joint_prob_mat(:,:,4) = max(joint_prob_mat(:,:,4),joint_prob_mat(:,:,4)');

                end % if adjacent vertices
            end
        end

        joint_prob_mat = (joint_prob_mat + TOL) ./ (1 + 4*TOL);
        % Now just get the mutual informations from the joint probability matrices:

        bnets_MI_mats(:,:,i) = - repmat(marginal_probs_vec .* log2(marginal_probs_vec) + ...
            (1-marginal_probs_vec) .* log2(1-marginal_probs_vec), n, 1);
        bnets_MI_mats(:,:,i) = bnets_MI_mats(:,:,i) + (1-eye(n)).*bnets_MI_mats(:,:,i)';
        bnets_MI_mats(:,:,i) = bnets_MI_mats(:,:,i) + joint_prob_mat(:,:,1) .* log2(joint_prob_mat(:,:,1)) + ...
            joint_prob_mat(:,:,2) .* log2(joint_prob_mat(:,:,2)) + ...
            joint_prob_mat(:,:,3) .* log2(joint_prob_mat(:,:,3)) + ...
            joint_prob_mat(:,:,4) .* log2(joint_prob_mat(:,:,4));

        bnets_MI_mats(:,:,i) = bnets_MI_mats(bnets{i}.rev_order,bnets{i}.rev_order,i); % Transfer ordering back
            
    end

end % way if



% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




