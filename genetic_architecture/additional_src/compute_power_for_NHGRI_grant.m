% Do some power calculations for NHGRI grant 
a_vec = 2.5 * 10.^(-6:-3); % false positive probability
b_vec = 0.1:0.1:0.9; % false negative probability

NCP = zeros(length(a_vec), length(b_vec)); 
for i=1:length(a_vec)
    for j=1:length(b_vec)
        NCP(i,j) = two_types_errors_to_non_centrality_parameter(a_vec(i), b_vec(j));
    end
end

R = [[0 b_vec]' [a_vec' NCP]'];

savecellfile(num2str_cell(num2cell(R), 3),   'n_a_b_table_for_NHGRI.txt', [], 1);
