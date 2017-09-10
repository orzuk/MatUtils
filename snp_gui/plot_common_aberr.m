function plot_common_aberr(del_stretch_cell, amp_stretch_cell, chr_num, loc_start, loc_end, ...
    genes_chr, genes_loc_start, gene_symbols, exon_start_cell, exon_end_cell, snp_chr_vec, snp_loc_vec)


% extract sub-stretch cell according to loc_start and loc_end

del_stretches_exist = 0;
amp_stretches_exist = 0;
[sub_stretch_cell_del, samples_del] = extract_sub_stretch_cell(del_stretch_cell, chr_num, loc_start, loc_end);
if(length(sub_stretch_cell_del))
    del_stretches_exist = 1;
    samples = samples_del;
end

[sub_stretch_cell_amp, samples_amp] = extract_sub_stretch_cell(amp_stretch_cell, chr_num, loc_start, loc_end);
if(length(sub_stretch_cell_amp))
    amp_stretches_exist = 1;
    if(~del_stretches_exist)
        samples = samples_amp;
    end
end
if(amp_stretches_exist | del_stretches_exist)
    num_samples = length(samples);
    samples = sort(samples);
    samples = samples(end:-1:1);

    % plot gene exons locations
    [fig_ind, xmin, xmax, ymin, ymax, legend_vec, legend_cell]= plot_genes_exons...
        (chr_num, loc_start, loc_end, genes_chr, gene_symbols, genes_loc_start, exon_start_cell, exon_end_cell);

    % plot snps locations
    idx=find(snp_chr_vec==chr_num & snp_loc_vec >= loc_start & snp_loc_vec <= loc_end);
    set(gca,'XTick',snp_loc_vec(idx), 'XTickLabel','', 'FontSize', 6);

    if(del_stretches_exist)
        [del_start_end_cell, del_loh_cell, ret_snp_loc_vec, big_aberr_vec1] = get_stretches_limits(sub_stretch_cell_del, samples);
        big_aberr_vec = big_aberr_vec1;
    end
    if(amp_stretches_exist)
        [amp_start_end_cell, amp_loh_cell, ret_snp_loc_vec, big_aberr_vec2] = get_stretches_limits(sub_stretch_cell_amp, samples);
    end
    if(~del_stretches_exist)
        big_aberr_vec = big_aberr_vec2;
    end

    stretch_legend_flag = 0;
    stretch_legend_flag2 = 0;
    loh_legend_flag = 0;

    % plot samples stretches
    samples_y_loc = [ymin+0.02:(ymax-0.02-(ymin+0.02))/(num_samples-1) :ymax-0.02];
    set(gca,'YTick',samples_y_loc);
    set(gca,'YTickLabel',samples, 'FontSize', 7);

    legend_ind = length(legend_vec)+1;

    for a = 1:num_samples
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First Plot DEL stretches
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(del_stretches_exist & length(del_start_end_cell)>0)
            stretch_color = 'b';
            start_vec = del_start_end_cell{a,1};
            end_vec = del_start_end_cell{a,2};
            stretch_len_vec = del_start_end_cell{a,3};
            num_stretches = length(start_vec);
            % plot big aberrations
            if(length(big_aberr_vec))
                sample_big_aberr = big_aberr_vec(a);
                if(sample_big_aberr~=0)
                    if(sample_big_aberr==1) big_aberr_str = 'AMP'; end
                    if(sample_big_aberr==-1) big_aberr_str = 'DEL'; end
                    plot([xmin xmax], [samples_y_loc(a) samples_y_loc(a)], '-k');
                    text(xmin+(xmax-xmin)/2, ...
                        samples_y_loc(a)+0.01, big_aberr_str, 'FontSize', 6);
                end
            end
            if(num_stretches>0)
                for b = 1:num_stretches
                    plot([start_vec(b) end_vec(b)], [samples_y_loc(a) samples_y_loc(a)], ['-' stretch_color]);
                    ind_start = plot([start_vec(b) start_vec(b)], [samples_y_loc(a) samples_y_loc(a)], ['.' stretch_color]);
                    plot([end_vec(b) end_vec(b)], [samples_y_loc(a) samples_y_loc(a)], ['.' stretch_color]);
                    text(start_vec(b), samples_y_loc(a)+0.015, [num2str(stretch_len_vec(b))], ...
                        'FontSize', 7);
                    if(~stretch_legend_flag2)
                        legend_vec(legend_ind) = ind_start;
                        legend_cell{legend_ind} = ['Del. Stretch Start-End'];
                        legend_ind = legend_ind+1;
                        stretch_legend_flag2 = 1;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now Plot AMP stretches
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(amp_stretches_exist & length(amp_start_end_cell)>0)
            stretch_color = 'r';
            start_vec = amp_start_end_cell{a,1};
            end_vec = amp_start_end_cell{a,2};
            stretch_len_vec = amp_start_end_cell{a,3};
            num_stretches = length(start_vec);
            if(num_stretches>0)
                for b = 1:num_stretches
                    plot([start_vec(b) end_vec(b)], [samples_y_loc(a) samples_y_loc(a)], ['-' stretch_color]);
                    ind_start = plot([start_vec(b) start_vec(b)], [samples_y_loc(a) samples_y_loc(a)], ['.' stretch_color]);
                    plot([end_vec(b) end_vec(b)], [samples_y_loc(a) samples_y_loc(a)], ['.' stretch_color]);
                    text(start_vec(b), samples_y_loc(a)+0.015, [num2str(stretch_len_vec(b)) ], ...
                        'FontSize', 7);
                    if(~stretch_legend_flag)
                        legend_vec(legend_ind) = ind_start;
                        legend_cell{legend_ind} = 'Amp Stretch Start-End';
                        legend_ind = legend_ind+1;
                        stretch_legend_flag = 1;
                    end
                end
            end
        end

        % plot LOH del
        if(del_stretches_exist & length(del_loh_cell)>0)
            sample_loh_vec = del_loh_cell{a,1};
            num_loh = length(sample_loh_vec);
            if(num_loh>0)
                for b = 1:num_loh
                    ind_start = plot([sample_loh_vec(b) sample_loh_vec(b)], [samples_y_loc(a) samples_y_loc(a)], '*m');
                    if(~loh_legend_flag)
                        legend_vec(legend_ind) = ind_start;
                        legend_cell{legend_ind} = 'LOH';
                        legend_ind = legend_ind+1;
                        loh_legend_flag = 1;
                    end
                end
            end
        end
        if(amp_stretches_exist & length(amp_loh_cell)>0)
            % plot LOH amp
            sample_loh_vec = amp_loh_cell{a,1};
            num_loh = length(sample_loh_vec);
            if(num_loh>0)
                for b = 1:num_loh
                    ind_start = plot([sample_loh_vec(b) sample_loh_vec(b)], [samples_y_loc(a) samples_y_loc(a)], '*m');
                    if(~loh_legend_flag)
                        legend_vec(legend_ind) = ind_start;
                        legend_cell{legend_ind} = 'LOH';
                        legend_ind = legend_ind+1;
                        loh_legend_flag = 1;
                    end
                end
            end
        end
    end


    title(['Stretches : Chr' num2str(chr_num) ':' num2str(loc_start) ' - ' num2str(loc_end)],'FontSize', 10);
    legend(legend_vec, legend_cell, 'Location', 'Best', 'FontSize', 7);
end