% Do a strupid check to see if the very very first integral is correct !
% We want to keep the top out of three

iters = 1000000;

c1_vec = [0:0.1:2];

prob_int = zeros(1, length(c1_vec));
prob_sim = zeros(1, length(c1_vec));
i=1;
for c1=c1_vec
    C_vec = [c1, 0, -1];


    %%% First do simulations
    noisy_C_vec = repmat(C_vec, iters, 1) + randn(iters, 3);

    [max_val max_ind ] = max(noisy_C_vec, [], 2);

    % output the maximal correlation
    prob_sim(i) = sum(max_ind == 1)/iters;

    cur_p_sim = prob_sim(i)

    %%% Now do the calculation using the stupid integral and see if it is the
    %%% same ...
    P_f1 = quadl('P_f1_integrand', -999, 999, [], [], C_vec)
    P_f0 = quadl('P_f0_integrand', -999, 999, [], [], C_vec);

    prob_integral =P_f1 / (P_f0 + P_f1)

    prob_int(i) = P_f1;
    i=i+1
end


figure; hold on; title('Compare simulated and integral'); xlabel('C1'); ylabel('Prob');
plot(c1_vec, prob_sim); plot(c1_vec, prob_int, 'r'); legend('sim', 'int');