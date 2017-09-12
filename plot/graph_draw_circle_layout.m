% Draw graph in a circular layout 
% 
% Input: 
% A - graph adjancancy matrix
% 
function graph_draw_circle_layout(A)

n = length(A); r = 0.4;
theta_vec = (1:n) .* (2*pi/n);
x_vec = 0.5 + r.*sin(theta_vec); y_vec = 0.5 + r.*cos(theta_vec); 

graph_draw_given_layout(A, [x_vec' y_vec']'); 

