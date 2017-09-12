% Simulate the Gladiators riddle
%
% Input: 
% n - initial number of gladiators
% iters - number of times to play the game
% 
% Output: 
% num_left - average number of gladiators left for the winning side 
% 
function [num_left gap_vec] = gladiators(n, iters)

num_left = 0; gap_vec = zeros(2*n,1); 
for j=1:iters
    if(mod(j, 10) == 0)
        run_iter = j
    end
    tmp_gap_vec = zeros(2*n,1);
    a=n; b=n; % start both sides with equal number 
    r = rand(1,2*n); % randomize all rounds in advance 

    for i=1:2*n % loop on fights 
        if(r(i) < a / (a+b))
            b=b-1;
        else
            a=a-1;
        end
        if(a == 0)
            num_left = num_left + b;
            break;
        end
        if(b == 0)
            num_left = num_left + a;
            break;
        end
        tmp_gap_vec(i) = abs(b-a); 
    end
    tmp_gap_vec(i:end) = abs(b-a);  % keep bias
    
    gap_vec = gap_vec + tmp_gap_vec; 
end
num_left = num_left / iters; % normalize
gap_vec = gap_vec ./ iters; 
