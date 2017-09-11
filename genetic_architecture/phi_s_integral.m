% Indefinite integral of phi_s times t^k dt
%
% Input:
% x - point in which to evaluate integral
% S - selection coefficient (in units of S=4Ns. Assume it is NEGATIVE!!!)
% phi_order - power of x
% grr_minus_one - additional parameter (optional) for disease
%
% Output:
% r - int_x t^k phi_s(t) dt. This is an indefinite integral!!
%
function r = phi_s_integral(x, S, phi_order, grr_minus_one)

if(~exist('phi_order', 'var'))
    phi_order = 0;
end
if(length(S) == 1)
    S = repmat(S, size(x));
end
if(length(x) == 1)
    x = repmat(x, size(S));
end
switch phi_order
    case 0 % int phi_s(t) dt
        r = ( Ei(S.*(x-1)) - log(1-x) + log(x) - exp(-S).*Ei(S.*x) ) ./ (1-exp(-S));
        r(S == 0) = log(x(S==0));
    case 1 % int t*phi_s(t) dt
        r = ( Ei(S.*(x-1)) - log(1-x) ) ./ (1-exp(-S));
        if(~isempty(find(S == 0,1)))
            r(S == 0) = x(S == 0);
        end
        null_inds = find(isnan(r)); % here S is too big - use asymptotics.
        if(~isempty(null_inds))   % Wrong Asymptotics !!!
            r(null_inds) = exp(S(null_inds) .* x(null_inds)) ./ (S(null_inds).*(1-x(null_inds))); % before we divided by S(null_inds)
        end
    case 2 % int t^2*phi_s(t) dt
        r = ( Ei(S.*(x-1)) - log(1-x) - x + exp(S.*(x-1))./S ) ./ (1-exp(-S));
        r(S == 0) = x(S == 0).^2./2;
        null_inds = find(isnan(r)); % here S is too big - use asymptotics.
        r(null_inds) = -exp(-S(null_inds) .* x(null_inds)) ./ S(null_inds).^2;
    case {-1, 'var', 'het'}  % int t*(1-t)*phi_s(t) dt (heterozygosity)
        r = (x - exp(-S.*(1-x))./S) ./ (1-exp(-S)); % Missing here !!!!
        if(~isempty(find(S == 0,1)))
            r(S == 0) = x(S == 0)-x(S == 0).^2./2;
        end        
    case -2 % int t^2*(1-t)*phi_s(t) dt  (this gives mean heterozygosity) 
        r = ( x.^2./2 - exp(S.*(x-1)) .* (x.*S-1) ./ S.^2 ) ./ (1-exp(-S));
        if(~isempty(find(S == 0,1)))
            r(S == 0) = x(S == 0).^2./2 - x(S == 0).^3 ./ 3;
        end
    case 'disease' % mean variance explained for disease: % int t*(1-t)*phi_s(t) dt 
        r = -exp(-S.*(grr_minus_one+1)./grr_minus_one) .* ...
            ( (grr_minus_one.*x+1) .* S .* real(Ei(S .* (x  + 1/grr_minus_one))) + ...
            grr_minus_one .* exp(S ./ grr_minus_one) .* (exp(S) - exp(S.*x)) ) ./ ...
            ((1-exp(-S)).* grr_minus_one.^2 .* (grr_minus_one .* x + 1));
%         r2 = -exp(-S.*(grr_minus_one+1)./grr_minus_one) .* ...
%             ( (grr_minus_one.*x+1) .* S .* real(exp(S .* (x  + 1/grr_minus_one)) ./ (S .* (x  + 1/grr_minus_one))) + ...
%             grr_minus_one .* exp(S ./ grr_minus_one) .* (exp(S) - exp(S.*x)) ) ./ ...
%             ((1-exp(-S)).* grr_minus_one.^2 .* (grr_minus_one .* x + 1));
%         r_num =  ( (grr_minus_one.*x+1) .* S .* real(Ei(S .* (x  + 1/grr_minus_one))) + ...
%             grr_minus_one .* exp(S ./ grr_minus_one) .* (exp(S) - exp(S.*x)) ) 
%         r2_num = ( (grr_minus_one.*x+1) .* S .* real(exp(S .* (x  + 1/grr_minus_one)) ./ (S .* (x  + 1/grr_minus_one))) + ...
%             grr_minus_one .* exp(S ./ grr_minus_one) .* (exp(S) - exp(S.*x)) ) 
%         ttt = -grr_minus_one ./ ...
%             ((1-exp(-S)).* grr_minus_one .* (grr_minus_one .* x + 1));
%         log(r(2:end)), log(r2(2:end)), log(ttt(2:end))
%         exact_ei = log(-real(Ei(S .* (x  + 1/grr_minus_one))))
%         approx_ei = log(-exp(S .* (x  + 1/grr_minus_one)) ./ (S .* (x  + 1/grr_minus_one)))
%         
%         log(r_num(2:end)), log(r2_num(2:end))
        
        % % %         (grr_minus_one .* exp(S.*(x-1)) - ...
% % %             S .* (grr_minus_one.*x-1) .* exp (S .* (1./grr_minus_one-1)) .* Ei(S .* (grr_minus_one.*x-1) ./ grr_minus_one)) ./ ...
% % %             ((1-exp(-S)).* grr_minus_one.^2 .* (grr_minus_one .* x + 1));
        if(~isempty(find(S == 0,1)))
            r(S == 0) = -((grr_minus_one+1) ./ (grr_minus_one .* x(S == 0) + 1) + log(grr_minus_one .* x(S == 0) + 1)) / grr_minus_one^2; % ??? Wrong !!!!
        end
        null_inds = find(isnan(r) | isinf(r)); % here S is too big - use asymptotics.
        if(~isempty(null_inds)) % use 2nd-order asymtotics of Ei
            r(null_inds) = -( 1 + grr_minus_one ./ (S(null_inds) .* (grr_minus_one .* x(null_inds) + 1)) .* ...
                exp(S(null_inds) .* (x(null_inds)-1)) ) ./ ...
                ((1-exp(-S(null_inds))) .* grr_minus_one .* (grr_minus_one .* x(null_inds) + 1));
        end
end % switch


%r = (1-exp(-S*(1-x))) / (x.*(1-x).*(1-exp(-S)));
