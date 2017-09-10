% Calculate the angle between two vectors 
function alpha = two_vecs_to_angle(V1, V2, deg_flag, varargin)


alpha = acos( max(min(V1*V2' ./ ( repmat(sqrt(sum(V1.^2,2)), 1, size(V2,1)) .* ...
repmat(sqrt(sum(V2.^2,2)), 1, size(V1,1))' ), 1), -1) );


% This is slower but maybe more accurate: 
% norm_v1_vec = zeros(size(V1, 1), 1); norm_v2_vec = zeros(size(V2, 1),1); 
% for i=1:size(V1, 1)
%     norm_v1_vec(i) = norm(V1(i,:));
% end
% for i=1:size(V2, 1)
%     norm_v2_vec(i) = norm(V2(i,:));
% end
% alpha = acos( max(min(V1*V2' ./ ( repmat(norm_v1_vec, 1, size(V2,1)) .* ...
% repmat(norm_v2_vec, 1, size(V1,1))' ), 1), -1) );



% correct from radians to degrees
if(nargin > 2)
    if(deg_flag == 1) % degrees
        alpha = alpha .* (180 / pi); %  57.2957795; % conversion from radians to degrees
    end
end