% Correct data of two expression vectors using lowess algorithm (by Amit Zeisel)
%
% Input: x, y expression vector befor log, there is no
% meaning for the order.
%
% Output: X, Y vectors after correction, R vector
% of the local correction.
%
% Written by Amit Zeisel
%
function [DATA,R]=lowess_two_arrays(data)
x=data(:,1);
y=data(:,2);
j=0;
jump=1;
a=log2(x.*y);%new variable
b=log2(x./y);%new variable
N=length(a);
[a_sort,XI]=sort(a);%transform to the sorted vectors
b_sort=b(XI);
in_same=find(a_sort==a_sort(1)) ;
if in_same>0
    a_sort=a_sort(max(in_same)+1:end);
    b_sort=b_sort(max(in_same)+1:end);
end
n=length(a_sort);
n_intervals=70;
dx=(a_sort(end)-a_sort(1))/n_intervals;%half of the width of the group for each point by values
% X=x(XI); Y=y(XI);
R=ones(floor(n/10),1);%intiate R
w=round(0.025*n);%half of the width of the group for each point by number of elements
dif=-a_sort(1:n-2*w)+a_sort(2*w:end-1);%the vector of differences if taking th width by number of elements
if length(find(dif<dx))<0.05*length(dif)%if more than 3% of the differences smaller than dx,
    %     the smoothing is one by the width of the group by values, if not, by
    %     number of elements
    for i=1:jump:n%make the calculation for one of ten elements
        j=j+1;
        clear B
        if i>w+1 & i<n-w
            R(j)=polyval(polyfit(a_sort(i-w:i+w),b_sort(i-w:i+w),1),a_sort(i));%mean(b_sort(i-w:i+w));%
        elseif i<=w
            R(j)=polyval(polyfit(a_sort(1:i+w),b_sort(1:i+w),1),a_sort(i));%mean(b_sort(1:i));%
        else
            R(j)=polyval(polyfit(a_sort(i-w:end),b_sort(i-w:end),1),a_sort(i));%mean(b_sort(i:end));%
        end
    end
else
    for i=1:jump:n
        j=j+1;
        clear B
        in=find(a_sort<a_sort(i)+dx & a_sort>a_sort(i)-dx);%find the elements in the in the group
        %         B=b_sort(in);
        R(j)=polyval(polyfit(a_sort(in),b_sort(in),1),a_sort(i));%mean(B);%
    end
end
% R=2*R;%*interp1(a_sort(1:jump:end),R,a_sort,'spline');%interpolate R
Y=2.^((a_sort-b_sort+R)/2);%calculate the correct variables
X=Y.*2.^(b_sort-R);
X=[x(XI(1:max(in_same)));X];
Y=[y(XI(1:max(in_same)));Y];
R=[ones(max(in_same),1);R];
%return to the original order of elements
X(XI)=X;
Y(XI)=Y;
R(XI)=R;
R=2.^R;
DATA=[X,Y];
