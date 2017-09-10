% Calculate correction factors for the IBH sum procedure
% according to description in the supp.info.
% It uses the Gaussian approximation and hence valid for m>100

function [C,s]=get_IBHsum_fac(m)

s1 = round(m*0.1);
s2 = round(m*0.35);
C = C_m_m(m,s1);
s = secant11(C,s1,s2,m);
% s = secant(s1,s2,m);

end

% Auxiliary function for ..
function f=C_m_m(m,s)
m0=m;
f=zeros(1,length(m0));
for i=1:length(m0)
    sigma = sqrt((m0(i)-1)/3);
    mult = 5;
    if (m0(i)-1 + mult*sigma <= s)
        integr = 0;
    else
        t =max(s, m0(i)-1-mult*sigma):(1e-3):min(m, m0(i)-1+mult*sigma);
        int_arg=normpdf(t,m0(i)-1,sqrt((m0(i)-1)/3))./t;
        integr = trapz(t,int_arg);
    end
    
    f(i)=m0(i)*(normcdf(s,m0(i)-1,sqrt((m0(i)-1)/3))/s + ...
        integr+(1-normcdf(m,m0(i)-1,sqrt((m0(i)-1)/3)))/m);
end
end

% Auxiliary function for ..
function f=C_m(m,m0,s)
f=zeros(1,length(m0));
for i=1:length(m0)
    sigma = sqrt((m0(i)-1)/3);
    mult = 5;
    if (m0(i)-1 + mult*sigma <= s)
        integr = 0;
    else
        t =max(s, m0(i)-1-mult*sigma):(sigma*1e-4):min(m, m0(i)-1+mult*sigma);
        int_arg=normpdf(t,m0(i)-1,sqrt((m0(i)-1)/3))./t;
        integr = trapz(t,int_arg);
    end
    
    f(i)=m0(i)*(normcdf(s,m0(i)-1,sqrt((m0(i)-1)/3))/s+integr+(1-normcdf(m,m0(i)-1,sqrt((m0(i)-1)/3)))/m);
end
end

% Auxiliary function for ..
function res = secant11(C,s1,s2,m)

x1=s1;
f1=m0maxC(m,x1)-C;
x2=s2;
f2=m0maxC(m,x2)-C;
d=abs(x2-x1);
while d>1
    d
    xnew=round(x2-(x2-x1)/(f2-f1)*f2);
    if xnew<=0
        xnew=max([x1-round(0.1*m),round(0.01*m)]);
    end
    x1=x2;
    x2=xnew;
    f1=m0maxC(m,x1)-C;
    f2=m0maxC(m,x2)-C;
    d=abs(x2-x1);
end
res=x2;
end


% Auxiliary function for ..
function [t,m0]=m0maxC(m,s)
m1=max([(s-3*round(m/100)),5]):round(m/100):max([(s+10*round(m/100)),round(m*0.6)]);
f=C_m(m,m1,s);
[temp,im]=max(f);
x1=m1(max([im-1,1]));
x2=m1(min([im+1,length(m1)]));
d=x2-x1;
while d>2
    m1=x1:max([round((x2-x1)/100),1]):x2;
    f=C_m(m,m1,s);
    [t,im]=max(f);
    x1=m1(max([im-1,1]));
    x2=m1(min([im+1,length(m1)]));
    d=x2-x1;
end
m0=(x2+x1)/2;
t=C_m(m,m0,s);
end

%
% function res = secant(s1,s2,m)
% epsilon=1e-3;
% d=10;
% while (abs(d)>epsilon)
%     d = (s2-s1)/(G(m,s2) - G(m,s1)) *G(m,s2);
%
%     s1=round(s2);
%     s2=round(s2-d);
%     if (s1==s2)
%         break;
%     end
% end
% res=s2;
% end
%
%
% function res = G(m,s)
%
% C_last = C_m(m,m,1);
% m1 = 5;
% m2 = round(m1+m/10);
%
% diff_1 = C_m(m,m1+1,s) - C_m(m,m1,s);
% while (diff_1 < 0)
%     m1=round(m1+2*(m/10));
%     diff_1 = C_m(m,m1+1,s) - C_m(m,m1,s);
% end
%
% diff_2 = C_m(m,m2+1,s) - C_m(m,m2,s);
% i=0;
% while(diff_2>0 & m2<=m)
%     m2 = round(m2+1*(m/10));
%     i = i+1;
%     diff_2 = C_m(m,m2+1,s) - C_m(m,m2,s);
% end
%
% while (abs(m2-m1)>1)
%
%     m_new = round((m1+m2)/2);
%
%     diff_new = C_m(m, m_new+1,s) - C_m(m,m_new,s);
%     diff_2 = C_m(m,m2+1,s) - C_m(m,m2,s);
%     diff_1 = C_m(m,m1+1,s) - C_m(m,m1,s);
%
%     if (diff_new*diff_1 < 0)
%         m2 = m_new;
%     end
%
%     if (diff_new*diff_2 <0)
%         m1 = m_new;
%     end
% end
%
%
%
% if (m1 ~=m2)
%     t = max(C_m(m,m1,s),C_m(m,m2,s));
% end
% % r1=C_m(m,m1,s)
% % r2=C_m(m,m2,s)
% % m1
% % m2
% res = abs(C_last - t);
%
% end
