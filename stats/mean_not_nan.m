%function vec = mean_not_nan(mat)
function vec = mean_not_nan(mat)

vec = zeros(1, size(mat,2));
vec(:) = nan;

for i = 1:size(mat,2)
    nan_ind = isnan(mat(:,i));
    not_nan_ind = find(nan_ind==0);
    if(length(not_nan_ind)>0)
        vec(i) = mean(mat(not_nan_ind,i));
    end
end