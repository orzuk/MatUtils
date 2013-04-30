% Perform refinement of the maximal sequences
%
function [Ref num_funcs] = Refinement(max_tab)

if(iscell(max_tab)) % work on normal fans
    Ref = max_tab{1};

    for i=2:length(max_tab)
        Ref = union(Ref, max_tab{i}, 'rows');
    end
    num_funcs = length(Ref);

else % work on 2d tables
    if(ndims(max_tab) == 3)
        [k,m,n] = size(max_tab);
        max_tab = reshape(max_tab, k, m*n);
    end

    [Ref I J] = unique(max_tab', 'rows'); % perform unique
    num_funcs = length(I); % num of different inference functions
    Ref = reshape(J, m,n);
    figure; imagesc(Ref); colorbar; title(['totally ' num2str(num_funcs) ' inf. functions']);
    xlabel('\epsilon'); ylabel('p'); 
end




