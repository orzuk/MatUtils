% Convert intervals sparse form to a density 
function [x_vec density_vec] = ...
    intervals_to_density(start_pos, end_pos, vals) % transfer singleton heights to [start end height] format

x_vec = min(start_pos):max(end_pos); 
density_vec = zeros(length(x_vec),1); 
for i=1:length(start_pos)
    density_vec( start_pos(i):end_pos(i)) = vals(i);
end

