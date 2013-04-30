% Compute distances between many pairs of 'random' points (distributions).
% Points can also be chosen from models. 
% We look for trends in distances distributions 
% 
% Input: 
% max_N - number of (binary) random variables 
% iters - number of points 
% rand_flag - how to choose points: 1 - choose randomly, 0 - choose from model 
% 
% Output: 
% KL_mat - matrix of pairwise Kullback-Leibler distances
% 
function KL_mat = KL_distances_rand_points(max_N, iters, rand_flag)

% The values of N we enumterate
KL_mat = zeros(max_N, iters);

% First check random points
switch rand_flag

    case 0 % both points random
        for N=2:max_N
            P = -log(rand(iters,2^N)); P = P ./ repmat(sum(P,2),1,2^N);
            Q = -log(rand(iters,2^N)); Q = Q ./ repmat(sum(Q,2),1,2^N);
            KL_mat(N,:) = sum(P' .* log2(P')) - sum(P' .* log2(Q'));  % Now compute KL (or entropy ??)
        end

    case 0.5 % one random, one from a BNT
        for N=2:max_N
            P = zeros(iters,2^N);
            ttt = cputime;
            B = create_bounded_degree_bnet(N, 1, iters);
            B_create_t = cputime - ttt;
            ttt = cputime;
            P = my_bnets_to_probs(B);
            % %             for j=1:iters
            % %                 P(j,:) = zuk_bnet_to_probs(B{j});
            % %             end
            B_to_prob_t = cputime - ttt;

            Q = -log(rand(iters,2^N)); Q = Q ./ repmat(sum(Q,2),1,2^N);
            KL_mat(N,:) = sum(P' .* log2(P')) - sum(P' .* log2(Q'));  % Now compute KL (or entropy ??)
        end


    case 1 % find closest Cij      
        KL_mat(2:end,:) = KL_mat(2:end,:) + 9999999999;   % check distance from closest Cij:
        for N=2:max_N
            P = -log(rand(iters,2^N)); P = P ./ repmat(sum(P,2),1,2^N);
            C = create_random_DAG(N,1);     % generate the click
            Cij = get_all_Cij_for_dag(C);   % generaet all the possible Cij's
            % Now enumerate them and find the smallest distance
            % (Try just one graph)
            for j=1:length(Cij)
                cur_KL = relative_entropy_to_bnet(P',Cij{j});
                KL_mat(N,:) = min(KL_mat(N,:),cur_KL);
            end
        end


    case 2 % find distance to a random graph with bounded in-degree
        KL_mat = [];
        KL_mat{1} = zeros(max_N, iters);  KL_mat{2} = zeros(max_N, iters);
        KL_mat{1}(2:end,:) = KL_mat{1}(2:end,:) + 9999999999;
        KL_mat{2}(2:end,:) = KL_mat{2}(2:end,:) + 9999999999;
        for N=2:max_N
            P = -log(rand(iters,2^N)); P = P ./ repmat(sum(P,2),1,2^N);
            %%%%            [G MI_mat KL_mat(N,:)] = random_Tree_to_point(P); % closest_DAG_to_point(P, 1);
            [G1 G2 MI_mat KL_mat{1}(N,:) KL_mat{2}(N,:)]= closest_and_random_Trees_to_point(P);


            % % %             dags = create_bounded_degree_DAG(N, 1, iters); % Create random DAGs
            % % %             for j=1:iters
            % % %                 KL_mat(N,j) = relative_entropy_to_bnet(P(j,:)', dags(:,:,j)); % closest_DAG_to_point(P, 1);
            % % %             end
        end


    case 2.5 % find closest graph with bounded in-degree
        KL_mat(2:end,:) = KL_mat(2:end,:) + 9999999999;
        for N=2:max_N
            P = -log(rand(iters,2^N)); P = P ./ repmat(sum(P,2),1,2^N);
            [G MI_mat KL_mat(N,:)] = closest_Tree_to_point(P); % closest_DAG_to_point(P, 1);
        end


    case 3 % Random graph with bounded in-degree, point is also in a low dimensional model
        KL_mat = cell(3,1);
        KL_mat{1} = zeros(max_N, iters);  KL_mat{2} = zeros(max_N, iters);
        KL_mat{1}(3:end,:) = KL_mat{1}(3:end,:) + 9999999999;
        KL_mat{2}(3:end,:) = KL_mat{2}(3:end,:) + 9999999999;

        for N=3:max_N
            ttt = cputime;
            B = create_bounded_degree_bnet(N, 1, iters);
            MI_mats = bnets_to_pairwise_MI(B,1); % No need to pass through P here
            create_time = cputime - ttt
            
            ttt = cputime;
            for j=1:iters
%                ttt = cputime;
                RandTree = uniform_spanning_tree(N,0); Tree = B{j}.dag + B{j}.dag';
                KL_mat{1}(N,j) = sum(sum(MI_mats(:,:,j) .* (Tree - RandTree)));  % Get the random KL difference
                KL_mat{3}(N,j) = sum(sum(MI_mats(:,:,j) .* Tree));  % Get the random KL difference, no rand tree
                [NextTree KL_mat{2}(N,j)] = NextBestTree(MI_mats(:,:,j), Tree);
                
            end
            Trees_time = cputime - ttt
            
        end
        KL_mat{1} = KL_mat{1} ./ 2; % We took each edge two times

        % % %     case 4 % find closest graph with bounded in-degree, but point is also in a low dimensional model
        % % %         KL_mat(3:end,:) = KL_mat(3:end,:) + 9999999999;
        % % %         for N=3:max_N
        % % %             N
        % % %             B = create_bounded_degree_bnet(N, 1, iters);
        % % %             MI_mats = bnets_to_pairwise_MI(B);
        % % %             for j=1:iters
        % % %                 %                [Tree MI_mat TmpKL] = closest_Tree_to_point(P(j,:)); %%%closest_DAG_to_point(P, 1);
        % % %                 %                [NextTree KL_mat(N,j)] = NextBestTree(MI_mat, Tree);
        % % %                 [NextTree KL_mat(N,j)] = NextBestTree(MI_mats(:,:,j), B{j}.dag + B{j}.dag');
        % % %             end
        % % %         end



end


