% Compute the GC content of a pwm
function gc = gc_for_pwm(pwm)

if(iscell(pwm))
    n = length(pwm);
    gc = zeros(n,1);
    for i=1:n
        gc(i) = gc_for_pwm(pwm{i});
    end
else
    gc = sum(pwm, 2);
    gc = sum(gc(2:3)) / sum(gc);
end
