% Plot the sigmoid curve along with the risk values for different genotypes
function plot_sigmoid_architecture(params_vec, architecture_formula, plot_flag)

a = params_vec{min_ind}(end-2);
b = params_vec{min_ind}(end-1);
c = params_vec{min_ind}(end);
q_z  = (1-sqrt(1-4*z_std(min_ind)^2))/2; % perform convulution with a binary variable
tmp_x_vec = [-10:0.01:10];
tmp_y_vec = q_z + (1-2*q_z) * 0.5 *(1+ tanh(a.*(tmp_x_vec-b))); % ./ (a+b.*exp(-c.*tmp_x_vec));

z = genetic_architecture(x_vec, 'sigmoid', ...
    params_vec{min_ind}, z_std(min_ind), 1); % generate outputs
if(plot_flag)
    figure; hold on; plot(tmp_x_vec, tmp_y_vec);
    
    plot(sum(x_vec .* repmat(vec2row(params_vec{min_ind}(1:N)), M, 1), 2), z, 'r*');
    title(['Best sigmoid architecture. N=' num2str(N) ', f=' num2str(f(1),precision) ...
        ', freq.=' num2str(good_architectures(arch_ind).mu,precision) ...
        ', h_{add}=' num2str(good_architectures(arch_ind).h_add,precision) ...
        ', h=' num2str(good_architectures(arch_ind).h,precision) ...
        '. z = 0.5 *(1 + tanh (' num2str(a,precision) ...
        '(x-' num2str(b,precision) '))']); %  e^{-' num2str(c,precision) 'x})']);
    legend('sigmoid', 'genotype-values');
    %                        xlabel('additive-risk-factor');
    xlabel(['Arch.: ' architecture_formula{min_ind}]);
    ylabel('prob. disease');
    snapnow; close all;
    
    figure; hold on; % plot cumulatie risk
    [z_sorted sort_perm] = sort(z);
    plot(cumsum(p_x_vec(sort_perm)), z_sorted); title('sigmoid cumulative risk');
    xlabel('frac. populaion'); ylabel('risk (prob. disease)');
    snapnow; close all;
end % if plot
