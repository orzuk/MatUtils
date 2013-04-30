% Divide matrix to blocks
% Assume that matrix is ordered already according to the blocks
%
% Input: 
% M - block matrix 
% 
%
% block_vec - cell array with indices, each cell representing a block 
% 
function block_vec = divide_matrix_to_blocks(M)

n = length(M); % assuem square matrix
block_ind_vec = zeros(n,1);

block_vec = cell(n,1);
block_ctr=1;
for i=1:n
    if(mod(block_ctr,10) == 0)
        run_block = block_ctr, run_i = i
    end
    if(~block_ind_vec(i)) % only if this index isn't already part of a block
        cur_block_inds = i; new_block_inds = find(M(i,:)); 
        while(~isempty(setdiff(new_block_inds, cur_block_inds)))
           cur_block_inds = new_block_inds;
            [I J] = find(M(cur_block_inds,:));             
            new_block_inds = union(J, cur_block_inds); 
           
%        add_to_block = find(M(i,:));
%        cur_block = union(cur_block, add_to_block);
        end
        block_ind_vec(cur_block_inds) = 1; 
        block_vec{block_ctr} = cur_block_inds;
        block_ctr = block_ctr+1;
    end
end
block_vec = block_vec(1:block_ctr-1); 

