%function cdf_f_name = get_cdf_f_name(chip)
function cdf_f_name = get_cdf_f_name(chip)

cdf_f_name = ['..\database\'];
if(strcmp(lower(chip), 'xba'))
    cdf_f_name = [cdf_f_name 'Mapping50K_Xba240.CDF'];
end

if(strcmp(lower(chip), 'hind'))
    cdf_f_name = [cdf_f_name 'Mapping50K_Hind240.CDF'];
end