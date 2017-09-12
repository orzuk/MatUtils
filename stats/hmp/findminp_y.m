% Find p minimizing Y probability according to HMP
p = -0.1; N = 6;
res = 0.01;

eps_vec = -0.5:0.0025:0.5-res;
p_vec = eps_vec*0+p;

P_Y = HMP_ProbY(eps_vec, p_vec, N);

[MIN_Y MIN_IND] = min(P_Y, [],  2);


pos_ind = find(MIN_Y >= 0);


POS_MIN_IND = MIN_IND(pos_ind)-1;

errors = length(POS_MIN_IND)-length(find(POS_MIN_IND == 22)) - length(find(POS_MIN_IND == 11))  


% for i=1:length(eps_vec)
%     sprintf('%lx\n', MIN_IND(i)-1)
% end