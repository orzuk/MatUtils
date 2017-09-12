% Counts nucleotides for a pwm
function counts = basecount_for_pwm(pwm)

counts = sum(pwm, 2); 

