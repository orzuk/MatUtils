%function stretch_len_vec = snp_large_num_in_loh(loh_vec, stretch_len_vec)
function stretch_len_vec = snp_large_num_in_loh(loh_vec, stretch_len_vec)

num_snps = length(stretch_len_vec);
large_num = 100;
% don't remove stretches with LOH
%loh_vec
loh_ind = find(loh_vec==1);
if(length(loh_ind)>0)
    for j = 1:length(loh_ind)
        stop = 0;
        ind = loh_ind(j);
        val = stretch_len_vec(ind);
        stretch_len_vec(ind) = large_num;
        % step down and assign large number
        while(~stop)
            ind = ind-1;
            if(ind==0)
                stop = 1;
            else
                if(stretch_len_vec(ind)==val)
                    stretch_len_vec(ind) = large_num;
                else
                    stop = 1;
                end
            end
        end

        stop = 0;
        ind = loh_ind(j);
        % step up and assign large number
        while(~stop)
            ind = ind+1;
            if(ind==(num_snps+1))
                stop = 1;
            else
                if(stretch_len_vec(ind)==val)
                    stretch_len_vec(ind) = large_num;
                else
                    stop = 1;
                end
            end
        end
    end
end