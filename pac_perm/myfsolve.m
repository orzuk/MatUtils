% Solve a non-linear equation for the special saddle point case 
function [a,avg_x,another_end] = myfsolve(mu, sigma, prior, alpha, one_side_flag, end_x)
another_end=0;
avg_x=0;
start_x=0;
a = find_c_alpha(mu, sigma, prior, alpha, one_side_flag, end_x);
%a1=a+100;
if (a<0)
    while(abs(a)>1.e-7)
%        a
%        a1=a;
        avg_x=(end_x+start_x)/2;
        a = find_c_alpha(mu, sigma, prior, alpha, one_side_flag, avg_x);
        if(a<0)
           end_x=avg_x;
        else
          start_x=avg_x; 
        end
    end
else % a>0
    another_end=1;
    'choose larger end_x'
end