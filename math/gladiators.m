% Simulate the Gladiators riddle
% The input: 
% n - number of gladiators
% iters - number of time to plat the game
% 
% The output: 
% num_left - average number of gladiators left 
% 
function num_left = Gladiators(n, iters)

num_left = 0;

for j=1:iters

    a=n; b=n;

    r = rand(1,2*n);


    for i=1:2*n
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

    end

end

num_left = num_left / iters;

