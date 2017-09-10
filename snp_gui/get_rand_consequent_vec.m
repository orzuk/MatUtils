%function num_consequent_vec = get_rand_consequent_vec(num_loc, num_events, num_start, max_num_consequent)
function num_consequent_vec = get_rand_consequent_vec(num_loc, num_events, num_start, max_num_consequent)
 
p = num_events/num_loc;
num_consequent_vec = [];

i=1;
cont = 1;
while (cont & i <= max_num_consequent)
%for i = 1:max_num_consequent
    num_loc_i = num_loc-(i-1);
    num_i_consequent_events = round(p^i*(1-p)^2*num_loc_i);
    num_consequent_vec = [num_consequent_vec repmat(i,1, num_i_consequent_events)];
    i = i+1;
    if(~ num_i_consequent_events)
        cont = 0;
    end
end
if(nargin >= 3)
    num_consequent_vec(find(num_consequent_vec < num_start)) = [];
end
