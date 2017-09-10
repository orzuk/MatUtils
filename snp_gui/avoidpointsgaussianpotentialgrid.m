% The function tries to find the 'best' curve, i.e. the curve which is as
% far as possible from the points given in the input.
function P = AvoidPointsGaussianPotentialGrid(grid_x, grid_y, data_x, data_y, sigma)

for i=1:length(grid_x)
    for j=1:length(grid_y)
        P(i,j) =  -sum( exp( - ((grid_x(i)-data_x).^2+(grid_y(j)-data_y).^2) ./(2*sigma^2) ) ); % here we have the points (x,y)
%%        a = grid_x(i); b = grid_y(j);
%%        P(i,j) = -sum( exp( -(data_y-a*data_x.^2-b.*data_x)./(2*sigma^2) ) ); % Here we have the coefficients (a,b)
    end
end

figure; imagesc(P); colorbar;

