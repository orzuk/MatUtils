% PLot SFS obtained from Wright Fisher MOdel - New function ! 
function [plot_x_vec, plot_y_vec, legend_vec] = ...
    FisherWrightPlotResults(freq_struct, absorption_struct, simulation_struct, D, ...
    s, mu, num_bins, init_str, iters, ...
    fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, all_s_output_file, mathematica_flag) % plot result


plot_x_vec = [];


