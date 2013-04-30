% Compute the prob. of getting a double at last move in backgammon (IBM riddle July 2009)
% (simulation should match recurrance relation)
function [p p_req P] = backgammon_riddle(n, iters, p_double, num_players)

x = rand(iters*n*num_players,1) < p_double; ctr=1; p=0; % generate many random dice
for i=1:iters
    k = repmat(n,num_players,1); % k(2) = 20*n;
    while(min(k) > 0)
        ctr=ctr+num_players; % update counter
        k = k - (1+1*x(ctr:ctr+num_players-1));
        %         for j=1:num_players
        %             ctr=ctr+1;
        %             if(x(ctr))
        %                 k(j)=k(j)-4;
        %             else
        %                 k(j)=k(j)-2;
        %             end
        %             if(k(j) <= 0)
        %                 break;
        %             end
        %         end
    end
    f = find(k <= 0, 1); p = p+x(ctr+f-1); %   p = p+x(ctr);
end
p = p / iters;

% Now compute analytically
P = zeros(n);
P(1,:) = p_double; P(:,1) = p_double; P(2,1) = p_double*(2-p_double);
for j=2:n
    P(2,j) = p_double + (1-p_double)*P(j,1);
end
for i=3:n
    P(i,2) = p_double*P(2,i-2) + (1-p_double)*P(2,i-1);
end
for i=3:n
    for j=3:n
        P(i,j) = p_double^2 * P(i-2,j-2) + (1-p_double)^2 * P(i-1,j-1) + ...
            p_double * (1-p_double) * P(i-1,j-2) + p_double * (1-p_double) * P(i-2,j-1);
    end
end
p_req = P(n,n); 
