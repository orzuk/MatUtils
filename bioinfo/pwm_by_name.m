% Get a pwm from a struct by it's name 
function pwm = pwm_by_name(pwms, name_str)

pwm = pwms{strmatch(name_str, pwms(:,1)),2};


