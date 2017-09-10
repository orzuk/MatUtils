% Get a set of pwms and get a list of all unique columns
% up to some error eps, in a new pwm
function [unique_pwm unique_inds] = unique_pwms_columns(pwms, epsilon, revcomp_flag)

num_pwms = length(pwms); 

L = zeros(num_pwms, 1); n = 0;
for i=1:num_pwms
   L(i) = size(pwms{i},2);  
end
L = [0 cumsum(L')]; n=L(end);
unique_pwm = zeros(4,n); 
for i=1:num_pwms
   unique_pwm(:,L(i)+1:L(i+1)) = pwms{i};
end

if(exist('revcomp_flag', 'var') || revcomp_flag)
    unique_pwm = [unique_pwm pwmrcomplement(unique_pwm)];
end

% [unique_pwm unique_inds] = unique(unique_pwm', 'rows'); 
[unique_pwm unique_inds] = unique_with_tol(unique_pwm', epsilon);
% figure; plot3(uu(1,:), uu(2,:), uu(3,:), '.')

    



