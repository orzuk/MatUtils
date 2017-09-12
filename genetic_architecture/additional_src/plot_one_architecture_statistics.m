% Plot various statistic for an attempt at finding one architecture with
% many parameters
function plot_one_architecture_statistics(V, v_genetic, v_additive_explained, ...
    p_z_x_marginal, architecture_formula)

%    figure; plot(sum(v_marginal,2), v_both, '.'); title('marginal variance vs. all variance explained');
%    xlabel('V(y|x_1) + V(y|x_2)'); ylabel('V(y|x_1,x_2)');
%        v_additive_explained = sum(repmat(V, 1, N) - v_marginal, 2); % amount of variance explained by additive effects
%        v_additive_explained2 = sum(repmat(V2, 1, N) - v_marginal2, 2); % amount of variance explained by additive effects
figure; hold on; plot(v_additive_explained ./ V, v_genetic ./ V, '.'); % plot ALL functions
%plot(v_additive_explained2 ./ V2, v_genetic2 ./ V2, 'r.'); % plot ALL functions. Additive explained is problematic when sampling
%  plot(N*(v_exact-v_marginal_exact)  ./ v_exact, (v_exact-v_environment_exact) ./ v_exact, 'g*');

for i=1:iters
    p_01(i) = p_z_x_marginal{i}(1,2);
    p_11(i) = p_z_x_marginal{i}(1,4);
end
figure; plot(p_01, p_11, '.'); title('different p_{ij}');
xlabel('p_{01}'); ylabel('p_{11}');

plot(v_additive_explained(valid_inds) ./ V(valid_inds), v_genetic(valid_inds) ./ V(valid_inds), 'rx');

min_val_x = min(v_additive_explained./V);
max_val_x = max(v_additive_explained./V);
min_val_y = min(v_genetic./V);
max_val_y = max(v_genetic./V);
min_val = min(min_val_x, min_val_y);
max_val = max(max_val_x, max_val_y);

plot(min_val-0.01:0.01:max_val+0.01,min_val-0.01:0.01:max_val+0.01,'k'); % what's this line ?
title_str = ['additive fraction vs. all fraction of variance explained, N=' ...
    num2str(N) ' f=' num2str(f(1)) '. Tried ' num2str(iters) ... % num2str(iters) ' iterations sampled '
    ' ' architecture_str ' architectures'];
title(title_str);
%            xlabel('additive fraction: h_{add} = N-(V(z|x_1) + ... + V(z|x_n))/V(z)'); ylabel('all fraction: h = 1-V(z|x)/V(z)');
xlabel(['Best Arch.: ' architecture_formula{min_ind}]);
legend('full enumeration',  'interesting-functions (\mu,h big enough)', 'y=x'); % 'sampling',
%             axis([min_val_x, max_val_x, min_val_y, max_val_y]);
snapnow; close all;
