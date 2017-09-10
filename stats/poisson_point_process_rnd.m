% Generate data from a Poisson point process in the cube [0,1]^n
%
% Input:
% lambda - rate of the poisson process
% dim - dimension of unit cube
% iters - how many sets of points to simulate
%
% Output:
% poiss_points - a matrix with 'dim' columns with coordinates of simulated points, or a cell
%                array of such matrices (when iters>1)
% npoints - number of points in each iteration
%
function [poiss_points npoints] = poisson_point_process_rnd(lambda, dim, iters)

npoints = poissrnd(lambda, iters, 1);   % number of points is distributed Poisson(lambda)
poiss_points = rand(sum(npoints), dim); % simulate points uniformly

if(iters>1) % seperate to cell array - one element for each iteration
    poiss_points = vec2cell(vec2row(poiss_points), cumsum(npoints));
end
