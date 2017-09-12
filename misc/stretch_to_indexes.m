% Spread the bits of x in the places which are ones in i.
% This works only up to 2^32 - need to fix that!!!
%
% Input:
% x - value. Only i lsb are important
% i - indicators of bits to spread x on
%
% Output:
% x_i - x values on the i indices
%
function x_i = stretch_to_indexes(x, i)


L = ceil(log2(i+1));
if(L > 50) % long values: use vpi
    i = vpi(i);
end
if(isa(i, 'vpi') || isa(x, 'vpi'))
    vpi_flag = 1;
    x_i=vpi(zeros(size(x,1),size(x,2))); % Copy input size
else
    vpi_flag = 0;
    x_i=zeros(size(x,1),size(x,2)); % Copy input size
end

new_ver_flag = 1;
if(new_ver_flag)
    if(isa(i, 'vpi'))
        i_bin = vpi2bin(i);
    else
        i_bin = dec2bin(i);
    end
    i_bin = find(i_bin(end:-1:1) == '1'); % where to put i
    for j=1:size(x,2)
        if(isa(x, 'vpi'))
            for k=1:size(x,1)
                tmp(k,:) = vpi2bin(x(k,j), length(i_bin)); % just get binary (put zeros at start)
            end
        else
            tmp = dec2bin(x(:,j), length(i_bin)); % just get binary (put zeros at start)
        end
        tmp = tmp(:,end:-1:1); tmp = tmp(:,1:length(i_bin));
        cur_x = repmat('0', size(x,1), max(i_bin));
        cur_x(:,i_bin) = tmp; % tmp(end:-1:1);
        if(vpi_flag)
            for k=1:size(x,1)
                x_i(k,j) = bin2vpi(cur_x(k,end:-1:1));
            end
        else
            x_i(:,j) = bin2dec(cur_x(:,end:-1:1));
        end
    end
else % use old (working for short integers) version
    
    num_on = hamming_weight(i);
    on_ind = zeros(1,num_on); % Find the 'on' indices
    k=1;
    for j=0:L
        if(mod(floor(i / 2^j), 2) == 1)
            on_ind(k) = j; k=k+1;
        end
    end
    
    for j=0:num_on-1
        x_i = x_i + bitshift(mod(floor(x ./ 2^j), 2), on_ind(j+1));
    end
end % if new version


