% A small function to test whether it is possible to get different values
% of inheritance and lods-ratio for two markers
function test_inheritance_and_lods(iters, p)

if(~exist('p', 'var') || isempty(p))
    p = rand(iters, 4); p = p ./ repmat(sum(p,2), 1, 4);
end
lods_ratio = p(:,4) .* (p(:,1) + p(:,3)) ./ (p(:,3) .* (p(:,2) + p(:,4))); % Pr(z=1|x=1) / Pr(z=1|x=0)


V = (p(:,3) + p(:,4)) .* (p(:,1) + p(:,2)); % var(z)
p_z_one_given_x_one = p(:,4) ./ (p(:,2) + p(:,4)); 
p_z_one_given_x_zero = p(:,3) ./ (p(:,1) + p(:,3)); 

v_env = (p(:,2) + p(:,4)) .* p_z_one_given_x_one .* (1-p_z_one_given_x_one) + ...
(p(:,1) + p(:,3)) .* p_z_one_given_x_zero .* (1-p_z_one_given_x_zero)

h = 1 - v_env ./ V; 

good_inds = intersect(find(h < 0.02), find(lods_ratio < 1.5)); 
good_inds = intersect(good_inds, find(lods_ratio > 1.0));

weird_inds = intersect(good_inds, find(lods_ratio > 2));
p_weird = p(weird_inds,:)
V_weird = V(weird_inds)
v_env_weird = v_env(weird_inds)
h_weird = h(weird_inds)
LODS_ratio_weird = lods_ratio(weird_inds)

figure; plot(h(good_inds), lods_ratio(good_inds), '.'); 
title('inheritance vs. lods ratio'); 
xlabel('h'); ylabel('lods ratio'); 
axis([0 0.02 1 1.5]); 
