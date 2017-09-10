% Convert a density to a sparse form of intervals 
function [start_pos end_pos vals] = ...
    density_to_intervals(x_vec, density_vec) % transfer singleton heights to [start end height] format

diff_density_vec = diff(density_vec); % detect change points 
change_points = find(diff_density_vec); 

start_pos = x_vec([1 vec2row(change_points)+1]); 
end_pos = x_vec([vec2row(change_points) length(density_vec)]); 
vals = density_vec([1 vec2row(change_points)+1]); 

