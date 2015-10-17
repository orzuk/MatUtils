% Get the first significant digits of a number
function x_msb = get_first_digits(x, direction, num_digits)


L = floor(log10(x));
if(~exist('direction', 'var') || isempty(direction))
    direction = 'round';
end
switch direction
    case {'down', 'floor'}
        x_msb = floor(x/10^L) * 10^L;
    case {'up', 'ceil'}
        x_msb = ceil(x/10^L) * 10^L;
    case 'round'
        
        x_msb = round(x/10^L) * 10^L;
end
% s = num2str(x);
% first_ind = find( (s ~= '0') & (s ~= '.'), 1);
%
% x_msb = str2num(s(first_ind:first_ind+num_digits-1)) * 10^(1-first_ind);

