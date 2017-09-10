%function D = plot_two_birds( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

AssignGeneralConstants;
num_turns = 1000;
strategy_vec = {'smart', 'stupid'};
for payoff_dist = {'unif', 'norm'} % loop on different distributions
    for i=1:length(strategy_vec)
        [alpha_vec{i} x_vec{i}] = siton_simulate_payoffs(num_turns, payoff_dist{1}, 0, strategy_vec{i});
        switch strategy_vec{i}
            case 'stupid' % run the stupid strategy for different number of turns
                for j=1:num_turns
                    [tmp_alpha tmp_x_vec] = siton_simulate_payoffs(j, payoff_dist{1}, 0, strategy_vec{i});
                    x_vec{i}(j) = tmp_x_vec(end);
                end
        end
    end
    
    figure; hold on;
    for i=1:length(strategy_vec)
        plot(x_vec{i}, [color_vec(i) '.']);
    end
    legend(strategy_vec, 4); xlabel('Num. Turns'); ylabel('Payoff');
    title(['Expected payoff for different strategy for ' payoff_dist{1} ' payoff quality']);
end


