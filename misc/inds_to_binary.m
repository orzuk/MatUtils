% Convert indices to binary
function bin_vec = inds_to_binary(inds_vec)

if(isa(inds_vec, 'vpi')) % special variable arithmatic integers
    bin_vec = double(vpi2bin(inds_vec-1)) - double('0');    
else
    bin_vec = double(dec2bin(inds_vec-1)) - double('0');
end

