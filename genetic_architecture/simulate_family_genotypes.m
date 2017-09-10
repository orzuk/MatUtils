% Simulate a genotype vector for each member of a family
%
% Input:
% max_generations - number of generations in family
% f_vec - minor allele frequencies. (Can also be heritabilities when genotypes are Gaussians!)
% iters - number of families to simulate
% compute_method - (not implemented yet) new! default is sampling but now we can also enumerate all possible vectors and give them weights
% genotypes_type_flag - binary (default) or Gaussians
%
% Output:
% family_tree - tree representing family members
% family_genotype_vec - vector of genotypes for each family member:
%    (i,j,k) is iter i, the j-th locus, the k-th individual
% family_phenotype_vec - vector of phenotypes for each family member:
%    (i,k) is iter i, the k-th individual (this output is optional)
%
function  [family_tree family_genotype_vec family_phenotype_vec] = ...
    simulate_family_genotypes(max_generations, f_vec, iters, ...
    compute_method, genotypes_type_flag) % simulate genotype vectors for entire family

if(~exist('genotypes_type_flag', 'var') || isempty(genotypes_type_flag))
    genotypes_type_flag = 'binary';
end
family_tree = generate_family_tree(max_generations);
N = length(f_vec); % number of loci
M = length(family_tree); % number of vertices
family_genotype_vec = zeros(iters,N,M); % genotype vec for all variables

cur_parents = get_graph_sources(family_tree);
for i=1:length(cur_parents) % assign values to all parents
    switch lower(genotypes_type_flag)
        case {'binary', 'discrete'}
            family_genotype_vec(:,:,cur_parents(i)) = initilize_x_vec_constants(N, 0, ...
                vec2row(f_vec), 'sampling', iters); % take a random set of paternal genotypes
        case 'gaussian' % here 'genotypes' are Gaussians with variance given by f_vec
            h_x = f_vec; % input 'f_vec' is actually the heritability of each genotype
            family_genotype_vec(:,:,cur_parents(i)) = randn(iters,N) .* ...
                repmat(sqrt(h_x.*0.5), iters, 1);
    end % switch genotypes type
end
for i=1:max_generations % loop on generations 
    cur_children = get_children(cur_parents, family_tree);
    for j=1:length(cur_children)
        cur_parents = get_parents(cur_children(j), family_tree);
        
        w = rand(iters,N) < 0.5; % choose maternal/paternal genotype
        family_genotype_vec(:,1:2:end,cur_children(j)) = ...
            family_genotype_vec(:,1:2:end,cur_parents(1)) .* w(:,1:2:end) + ...
            family_genotype_vec(:,2:2:end,cur_parents(1)) .* (1-w(:,1:2:end)); % get maternal genotype
        %        w = rand(iters,N) < 0.5; % choose paternal genotype (same randomization)
        family_genotype_vec(:,2:2:end,cur_children(j)) = ...
            family_genotype_vec(:,1:2:end,cur_parents(2)) .* w(:,2:2:end) + ...
            family_genotype_vec(:,2:2:end,cur_parents(2)) .* (1-w(:,2:2:end)); % get paternal genotype
    end % loop on children
    cur_parents = cur_children;
end

if(nargout>2) % simulate also phenotye. Currently implement LP model     
    family_phenotype_vec = zeros(iters,M); 
    for i=1:M % loop on family member 
        family_phenotype_vec(:,i) = ... %             genetic_architecture( ... )
            max(reshape(family_genotype_vec(:,1:2:end,i) + ...
            family_genotype_vec(:,2:2:end,i),iters,N/2) + ...
            randn(iters,N/2) .* repmat(sqrt(1-h_x(1:2:end)), iters, 1),[],2); 
    end
end    
    