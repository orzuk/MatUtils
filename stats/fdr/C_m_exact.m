% Calculate correction factors C for the IBH sum procedure 
% according to description in the supp.info. 
% It uses the exact uniform-sum density function (no approximation) and hence 
% the calculation is very slow, not recommanded for m>15;
%
% Input: 
% m - number of procedures
% s - IBH multiplicative correction factor (s)
% 
% Output: 
% f - IBH additive correction factor (C)
% 
function C = C_m_exact(m,s)

m0=m;

t1=s:1e-3:m;
t2 = 0:1e-3:s;
t3 = m:1e-3:2*m;
int_arg1 = zeros(1, length(t1));
int_arg2 = zeros(1,length(t2));
int_arg3 = zeros(1,length(t3));
for i = 1:length(t1)
    int_arg1(i) = uniform_sum_pdf(t1(i),m0)/t1(i);
    
end
for i=1:length(t2)
    int_arg2(i) = uniform_sum_pdf(t2(i),m0);
end
for i=1:length(t3)
    int_arg3(i) = uniform_sum_pdf(t3(i),m0);
end

C = m0*(trapz(t2, int_arg2)/s + trapz(t1,int_arg1) + trapz(t3, int_arg3)/m);

