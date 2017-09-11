% New: alternative simulation of a population.
% How do we compute this?
%
% Input:
% num_founders - how many different founders for the population. Average IBD sharing is 1/num_founders
% num_generations - determine how many recombination break-points we expect (num_generations for each 'genome/chromosome')
% num_people - number of  individuals to simulate
% f_vec - allele frequencies of SNPs
% input_IBD_mat - (optional) take the IBD mat from input - simulate only SNPs 
% input_founder_identities - (optional) 
% input_breakpoints_vec - (optional) 
% input_SNP_founder_identity_mat - (optional) 
%
% Output:
% IBD_mat - matrix of IBD sharing between each pair of individuals
% estimated_IBD_mat - matrix of estimaetd IBD sharing based on simulated SNPs 
% SNP_mat - SNP matrix (generated according to the IBD-sharing matrix)
% founder_SNP_mat - SNP matrix for the founders (unrelated)
% SNP_founder_identity_mat - identity of founder for each individual and each SNP
% founder_identities - cell array with founder identities for each individual 
% breakpoints_vec - breakpoints of recombination 
%
function [IBD_mat estimated_IBD_mat SNP_mat ...
    founder_SNP_mat SNP_founder_identity_mat founder_identities breakpoints_vec] = ...
    simulate_IBD_blocks_sharing(num_founders, num_generations, num_people, f_vec, ...
    input_IBD_mat, input_founder_identities, input_breakpoints_vec, input_SNP_founder_identity_mat)

AssignGeneralConstants;
if(num_people <= 50) % plot only when simulating a few individuals 
    plot_flag = 1;
else
    plot_flag = 0;
end

num_snps = length(f_vec);
SNP_mat = zeros(num_snps, num_people);
%founder_SNP_mat = zeros(num_snps, num_founders);
founder_SNP_mat = rand(num_snps, num_founders) < repmat(f_vec, 1, num_founders); % Generate SNP matrix


if(exist('input_IBD_mat', 'var') && (~isempty(input_IBD_mat))) % Take IBD matrix from input
    IBD_mat = input_IBD_mat;
    founder_identities = input_founder_identities;
    breakpoints_vec = input_breakpoints_vec;
    SNP_founder_identity_mat = input_SNP_founder_identity_mat;
else % Compute IBD matrix 
    SNP_founder_identity_mat = SNP_mat;
    SNP_pos = (1:num_snps) ./ num_snps;
    genome_len = 10^8;
    founders_chr_vec = mat2vec(repmat(1:num_founders, num_snps, 1));
    founders_start_vec = round(mat2vec(repmat(SNP_pos, num_founders, 1)') * genome_len);
    founders_end_vec = founders_start_vec;

    
    lambda = num_generations; % Assume one recombination per generation 
    [breakpoints_vec npoints] = poisson_point_process_rnd(lambda, 1, num_people); % simulate break-points
    founder_identities = vec2cell(ceil(rand(1, sum(npoints+1)).*num_founders), cumsum(npoints+1)); % simualte colors    
    for i=1:num_people % sort points to enable faster sorting later
        breakpoints_vec{i} = sort(breakpoints_vec{i});
    end    
    IBD_mat = compute_IBD_mat_internal(breakpoints_vec, founder_identities, npoints, num_founders); % compute IBD sharing

    for i=1:num_people % loop on people to simulate SNPs
        if(mod(i, 50) == 0)
            simulate_SNPs_person = i
        end
        cur_person_chr_vec = founder_identities{i};
        cur_person_start_vec = sort(round([0 breakpoints_vec{i}] * genome_len));
        cur_person_end_vec = sort(round([breakpoints_vec{i} 1]  * genome_len));
        
        [intersect_chr_vec intersect_start_vec intersect_end_vec ...
            intersect_inds1 intersect_inds2 ...
            regions_len1 regions_len2 regions_inter_len] = ...
            intersect_genomic_regions(founders_chr_vec, founders_start_vec, founders_end_vec, ...
            cur_person_chr_vec, cur_person_start_vec, cur_person_end_vec);
        [intersect_start_vec sort_perm] = sort(intersect_start_vec);
        intersect_end_vec = intersect_end_vec(sort_perm);
        intersect_chr_vec = intersect_chr_vec(sort_perm);
        if(length(intersect_chr_vec) ~= num_snps)
            errr = 1322
            if(length(intersect_chr_vec) > num_snps)
                intersect_chr_vec = intersect_chr_vec(1:num_snps);
            else
                intersect_chr_vec = ...
                    [intersect_chr_vec repmat(1, 1, num_snps-length(intersect_chr_vec))];
            end
        end
        SNP_founder_identity_mat(:,i) = intersect_chr_vec; % determine for each snp and each person who's the founder
%        SNP_mat(:,i) = founder_SNP_mat(sub2ind([num_snps num_founders], ...
%            1:num_snps, min(intersect_chr_vec, num_founders)));
    end
end % simulate IBD mat

%SNP_mat2 = zeros(size(SNP_mat)); 
SNP_founder_identity_mat = min(SNP_founder_identity_mat, num_founders); 
for i=1:num_people % A different way - just copy the SNPs from the founders
    SNP_mat(:,i) = founder_SNP_mat(sub2ind([num_snps num_founders], ...
        1:num_snps, vec2row(SNP_founder_identity_mat(:,i))));
end    
    
if(exist('input_IBD_mat', 'var') && (~isempty(input_IBD_mat))) % estimate IBD mat only first time
    estimated_IBD_mat = genetic_relationship_matrix(SNP_mat, 'binary', f_vec);
else % no need to reconstruct
    estimated_IBD_mat = [];
end
if(plot_flag)
    figure; subplot(2,1,1); imagesc(SNP_founder_identity_mat'); colorbar; title('Founder Inheritance Regions');
    subplot(2,1,2);  imagesc(SNP_mat'); colorbar; title('Population SNPs');
    figure; imagesc(IBD_mat-eye(num_people)); colorbar; title('IBD-sharing matrix');
    if(num_people <=20)
        figure; hold on;
        for i=1:num_people
            intervals_plot(sort([0 breakpoints_vec{i}]), sort([breakpoints_vec{i} 1]), ...
                founder_identities{i} + i*0.1/num_people, color_vec(i));
        end
        legend(num2str((1:num_people)'));
    end
    %    estimated_IBD_mat = estimate_IBD_sharing_from_snps(SNP_mat, f_vec, 1);
    title('Estimated IBD-sharing matrix');
    IBD_mat = IBD_mat - eye(num_people);
    if(~isempty(estimated_IBD_mat))
        estimated_IBD_mat = estimated_IBD_mat - eye(num_people);
        figure; hold on; plot(IBD_mat(:), estimated_IBD_mat(:), '.');
        [C C_pval] = corr(IBD_mat(:), estimated_IBD_mat(:));
        beta = regress(estimated_IBD_mat(:), [ones(length(IBD_mat(:)), 1) IBD_mat(:)]);
        x_vec = 0:0.01:1;
        plot(x_vec, x_vec .* beta(2) + beta(1), 'r');
        title(['IBD estimation corr=' num2str(C,3) ', n_{SNPs}=' num2str(num_snps)]);
        xlabel('True IBD'); ylabel('Estimated IBD');
        figure;  imagesc(estimated_IBD_mat); colorbar;
    end
end % if plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal function: Compute matrix of IBD sharing from set of intervals for each individual
function IBD_mat = compute_IBD_mat_internal(breakpoints_vec, founder_identities, npoints, num_founders)

num_people = length(breakpoints_vec);
IBD_mat = zeros(num_people); % IBD_mat2 = IBD_mat;
for old_way=0:0
    time_way = cputime; 
    if(old_way)
        for i=1:num_people % Compute IBD overlap between each pair of individuals. Currently the heaviest part!!!!
            if(mod(i, 50) == 0)
                simulate_IBD_person = i
            end
            for j=i+1:num_people % loop on second person
                %        [cur_breakpoints I J] = union(breakpoints_vec{i}, breakpoints_vec{j}); % find union
                [cur_breakpoints sort_perm] = sort([breakpoints_vec{i} breakpoints_vec{j}]);
                inv_sort_perm = inv_perm(sort_perm);
                cur_breakpoints = [0 cur_breakpoints 1];
                cur_interval_lens = diff(cur_breakpoints);
                temp_identities1 = zeros(1, npoints(i)+npoints(j)+1); temp_identities2 = temp_identities1;
                I = sort(inv_sort_perm(1:npoints(i))); J = sort(inv_sort_perm(npoints(i)+1:end));
                I = [0 I]; J = [0 J];
                temp_identities1(I+1) =  founder_identities{i}; % (2:end);
                temp_identities2(J+1) =  founder_identities{j}; % (2:end);
                I = [I npoints(i)+npoints(j)+1];
                J = [J npoints(i)+npoints(j)+1];
                for k=1:npoints(i)+1 % loop on all break-points
                    temp_identities1(I(k)+1:I(k+1)) = temp_identities1(I(k)+1);
                end
                for k=1:npoints(j)+1  % loop on all break-points in second person
                    temp_identities2(J(k)+1:J(k+1)) = temp_identities2(J(k)+1);
                end
                IBD_mat(i,j) = sum(cur_interval_lens(temp_identities1 == temp_identities2));
            end
        end
        IBD_mat = IBD_mat + IBD_mat' + eye(num_people);
    else % new way: refinement. Should be faster.
        all_breakpoints = sort(unique(cell2vec(breakpoints_vec))); % generate a refinment of all points 
        all_intervals_lens = diff([0 all_breakpoints 1]);
        num_all_breakpoints = length(all_breakpoints);     
        IBD_by_founder_mat = zeros(num_people, num_all_breakpoints+1);
        for i=1:num_people % intersect refinement with each individual
            [inter_starts inter_ends I] = ...
                intervals_intersect([0 breakpoints_vec{i}]+eps, [breakpoints_vec{i} 1]-eps, ...
                [0 all_breakpoints], [all_breakpoints 1], [], 1);
            IBD_by_founder_mat(i,:) = founder_identities{i}(I); % we've found the founder
        end
        %            IBD_mat2 = IBD_mat + slmetric_pw(IBD_by_founder_mat, IBD_by_founder_mat, 'hamming', Q);
    
        computed_IBD_by_founder_mat = 999    
        
        for i=1:num_all_breakpoints+1 % loop on intervals
            for j=1:num_founders % loop on founders
                J = find(IBD_by_founder_mat(:,i) == j); 
                IBD_mat(J,J) = IBD_mat(J,J) + all_intervals_lens(i);
            end 
        end
        IBD_mat = IBD_mat - diag(diag(IBD_mat)) + eye(num_people);
    end % if new way to compute 
    time_way = cputime-time_way;
    sprintf(['Method ' num2str(old_way) ' time=' num2str(time_way) ' sec.'])
end


