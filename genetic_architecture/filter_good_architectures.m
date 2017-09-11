% Pass only architectures which satisfy according to certain constraints
% 
% Input: 
% candidate_architectures, ...
% h_interval - possible interval for heretability
% h_add_interval - possible interval for additive heretability 
% freq_interval - possible interval for disease frequency in population
% ratio_interval - possible interval for ratio h_add / h
% lods_interval  - possible interval for lods-ratio for marginal risks
% penetrance_interval - possible interval for penetrance (twin risk)
% architecture_str - take only architectures of a certain type 
% plot_flag - whether to plot stuff (what?)
%
% Output: 
% valid_inds - indices of good architectures 
% valid_inds_mat - ??? 
% interesting_flag - found at least one interesting architecture
% 
function  [valid_inds valid_inds_mat interesting_flag] = ...
    filter_good_architectures(candidate_architectures, ...
    h_interval, h_add_interval, freq_interval, ratio_interval, ...
    lods_interval, penetrance_interval, architecture_str, ...
    plot_flag)

if(~exist('plot_flag', 'var') || isempty(plot_flag))
    plot_flag = 0;
end
valid_inds_str = {'h', 'freq.', 'r', 'penet.', 'h_{add}', 'h_{ij}', 'lods'};

if(isfield(candidate_architectures(1), 'mu')) % candidate architectures
    valid_inds_vec{1} = intersect(find(candidate_architectures(1).h > h_interval(1)), ...
        find(candidate_architectures(1).h < h_interval(2)));
    valid_inds_vec{2} = intersect(find(candidate_architectures(1).mu > freq_interval(1)), ...
        find(candidate_architectures(1).mu < freq_interval(2)));
    valid_inds_vec{3} = intersect(find(candidate_architectures(1).v_additive_explained ./ ...
        candidate_architectures(1).v_genetic > ratio_interval(1)), ...
        find(candidate_architectures(1).v_additive_explained ./ ...
        candidate_architectures(1).v_genetic < ratio_interval(2)));
    valid_inds_vec{4} = intersect(find(candidate_architectures(1).penetrance > ...
        penetrance_interval(1)), ...
        find(candidate_architectures(1).penetrance < penetrance_interval(2))); % take only interesting cases
    valid_inds_vec{5} = intersect(find(candidate_architectures(1).h_add > h_add_interval(1)), ...
        find(candidate_architectures(1).h_add < h_add_interval(2))); % take only interesting cases
%    for i=1:length(candidate_architectures)
        max_pairwise_vec = max(candidate_architectures(length(candidate_architectures)).v_pairwise_explained(:));
%    end
    valid_inds_vec{6} = find(candidate_architectures(1).v_additive_explained > ...
        max_pairwise_vec);
    valid_inds_vec{7} = intersect(find(max(candidate_architectures(1).lods_ratio_marginal) > ...
        lods_interval(1)), ...
        find(max(candidate_architectures(1).lods_ratio_marginal) < lods_interval(2)));
else % good architectures
    valid_inds_vec = cell(7,1); 
    for i=1:length(candidate_architectures)
        valid_inds_vec{1}(i) = ~isempty(intersect(find(candidate_architectures(i).h > h_interval(1)), ...
            find(candidate_architectures(i).h < h_interval(2))));
        valid_inds_vec{2}(i) = ~isempty(intersect(find(candidate_architectures(i).freq > freq_interval(1)), ...
            find(candidate_architectures(i).freq < freq_interval(2))));
        valid_inds_vec{3}(i) = ~isempty(intersect(find(candidate_architectures(i).v_additive_explained ./ ...
            candidate_architectures(i).v_genetic > ratio_interval(1)), ...
            find(candidate_architectures(i).v_additive_explained ./ ...
            candidate_architectures(i).v_genetic < ratio_interval(2))));
        valid_inds_vec{4}(i) = ~isempty(intersect(find(candidate_architectures(i).penet > ...
            penetrance_interval(1)), ...
            find(candidate_architectures(i).penet < penetrance_interval(2)))); % take only interesting cases
        valid_inds_vec{5}(i) = ~isempty(intersect(find(candidate_architectures(i).h_add > h_add_interval(1)), ...
            find(candidate_architectures(i).h_add < h_add_interval(2)))); % take only interesting cases
%        for j=1length(candidate_architectures)
            max_pairwise_vec = max(candidate_architectures(i).v_pairwise_explained(:));
%        end
        valid_inds_vec{6}(i) = ~isempty(find(candidate_architectures(i).v_additive_explained > ...
            max_pairwise_vec));
        valid_inds_vec{7}(i) = ~isempty(intersect(find(max(candidate_architectures(i).L_i_full) > ...
            lods_interval(1)), ...
            find(max(candidate_architectures(i).L_i_full) < lods_interval(2))));
    end
    for i=1:length(valid_inds_vec)
        valid_inds_vec{i} = find(valid_inds_vec{i});
    end
end
num_valid_inds = length(valid_inds_vec);
if(exist('architecture_str', 'var') && (~isempty(architecture_str)))
    valid_arch_inds = strmatch(architecture_str, candidate_architectures(1).architecture_str);
    for i=1:num_valid_inds
        valid_inds_vec{i} = intersect(valid_inds_vec{i}, valid_arch_inds);
    end
end
valid_inds = valid_inds_vec{1};
for i=1:num_valid_inds
    valid_inds = intersect(valid_inds, valid_inds_vec{i});
end
valid_inds_mat = zeros(num_valid_inds);
for i=1:num_valid_inds
    for ii=1:num_valid_inds
        valid_inds_mat(i,ii) = length(intersect(valid_inds_vec{i}, valid_inds_vec{ii}));
    end
end
if(plot_flag) % plot how many arch. satisfy pairwise interval conditions 
    figure; imagesc_with_labels(valid_inds_mat, valid_inds_str, valid_inds_str);
    for i=1:num_valid_inds
        for ii=1:num_valid_inds
            text(i,num_valid_inds-ii+1,num2str(valid_inds_mat(i,ii)));
        end
    end
    title('Different conditions intersection');
    figure; plot(v_additive_explained ./ (N.*V), max(lods_ratio_marginal,[],2), 'r.');
    xlabel('h_{add}'); ylabel('lods-ratio');
end

if(isempty(valid_inds)) % didn't found any one matching the criteria
    %                valid_inds = find(v_genetic ./ V > min_h); % make at least heretability large
    if(isfield(candidate_architectures(1), 'mu')) % candidate architectures
        [~, valid_inds] = max(vec2column(candidate_architectures(1).mu) .* ...
            vec2column(candidate_architectures(1).h));
    else
        for i=1:length(candidate_architectures)
            val_to_max(i) = candidate_architectures(i).freq .* ...
                candidate_architectures(i).h;
        end
        [~, valid_inds] = max(val_to_max);
    end
    interesting_flag = 0;
else
    interesting_flag = 1;
    found_good_interaction = 1111
end

