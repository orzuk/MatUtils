% Simulate a tennis set / tie breaker
function p_win = tennis_set(n, p, q, iters)


r = rand(iters,n*10); 

tie_breaker = 1; 
if(~tie_breaker)
    r(:,1:2:end) = r(:,1:2:end) < p;
    r(:,2:2:end) = r(:,2:2:end) > q;
else
    r(:,[1:4:end 4:4:end]) = r(:,[1:4:end 4:4:end]) < p;
    r(:,[2:4:end 3:4:end]) = r(:,[2:4:end 3:4:end]) > q; 
end
x = cumsum(r,2); y = cumsum(1-r,2); 

winner = zeros(iters,1); % x1 = zeros(iters,1); x2 = zeros(iters,1);

for i=n:n*10
    u_inds = find(winner==0);
    f = intersect(u_inds,find(x(:,i)>=max(n,y(:,i)+2))); 
    g=intersect(u_inds,find(y(:,i)>=max(n,x(:,i)+2))); 
    winner(f) = 1; winner(g) = -1;
end

undecided_games = find(winner == 0)
p_win = mean((winner+1)./2); 


