% Calculates the distance of every permutation from the identity one, 
% by how many transpositions are needed
% 
% Input: 
% N - perumtations size
% 
function swaps_vector  = calc_swaps_distance(N)

num_swaps = N*(N-1)/2;

all_perms = perms(1:N);

swaps_vector = zeros(1,factorial(N)) + N^2; % give a ridiculus high number to begin with

cur_perms =1:N; % Take the identity permutation
for swaps=1:N-1 % We know that the maximal swaps number is N-1
    new_perms = zeros(size(cur_perms,1)*num_swaps,N);
    swap_ind=0;
    for i=1:N
        for j=i+1:N
            cur_swap_perm = 1:N; cur_swap_perm(i)=j; cur_swap_perm(j)=i;
%             size( cur_perms(:,cur_swap_perm))
%              (N*(i-1)+j)*size(cur_perms,1)+1
%              (N*(i-1)+j+1)*size(cur_perms,1)
%             size( new_perms(:, (N*(i-1)+j)*size(cur_perms,1)+1: (N*(i-1)+j+1)*size(cur_perms,1)) )
%             size(cur_perms,1)
%            swap_ind*size(cur_perms,1)+1
%             (swap_ind+1)*size(cur_perms,1)
%             size(cur_perms(:,cur_swap_perm))
%             size(new_perms(:, (swap_ind*size(cur_perms,1)+1: (swap_ind+1)*size(cur_perms,1))))
%             size(new_perms)
%             size(cur_perms)
            new_perms( (swap_ind*size(cur_perms,1)+1: (swap_ind+1)*size(cur_perms,1)),:) = cur_perms(:,cur_swap_perm);
            swap_ind=swap_ind+1;
        end
    end
    new_perms = unique(new_perms,'rows');
    cur_perms = new_perms;
  
%     size(all_perms)
%     size(new_perms)
    [vals new_perms_indexes] = intersect(all_perms,new_perms,'rows') ;
    
    swaps_vector(new_perms_indexes) = min(swaps,  swaps_vector(new_perms_indexes));

end
[vals unit_index] = intersect(all_perms, 1:N, 'rows');

swaps_vector(unit_index)=0; % Correct the unit permutation