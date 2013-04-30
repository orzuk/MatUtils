function [h h_lines min_intensity_seg_all max_intensity_seg_all]=drawChromMap(handles, is_empty)

min_intensity_seg_all=[];
max_intensity_seg_all=[];

%chr = handles.chr;
chr_num = handles.chr_num;

includeSexChr = 1; %length(handles.ChromsData)>22; %ismember(23,chr);

left=0.05;
bottom=0.25;
width=0.9;
height=0.73;

chr_change = find(diff(chr_num) ~= 0);
if ~includeSexChr
    chr_change = chr_change(1:21);
end

if ~is_empty
    
    if get(handles.checkbox_all_samples','value')
        not_empty_chrom=find(~cellfun('isempty',handles.ChromsData));

        %calculate min and max of intensity segments for all chromosomes
        min_intensity_seg_all=min(handles.ChromsData{not_empty_chrom(1)}.Segments(:,3));
        max_intensity_seg_all=max(handles.ChromsData{not_empty_chrom(1)}.Segments(:,3));

        for i=2:length(not_empty_chrom)
            min_intensity_seg_all=min([min_intensity_seg_all; handles.ChromsData{not_empty_chrom(i)}.Segments(:,3)]);
            max_intensity_seg_all=max([max_intensity_seg_all; handles.ChromsData{not_empty_chrom(i)}.Segments(:,3)]);
        end
    else
        min_intensity_seg_all=0; 
        max_intensity_seg_all=5; 
    end
    handles.min_intensity_seg_all=min_intensity_seg_all;
    handles.max_intensity_seg_all=max_intensity_seg_all;
end

if ~is_empty
%     axes(handles.axes_colorbar);
%     set(gca,'visible','on');
    colorbar_h=colorbar('peer',handles.axes_colorbar);
    colormap('default');
    ylim1=get(colorbar_h,'ylim');
    set(colorbar_h,'ytick',[ylim1(1) mean(ylim1) ylim1(2)]);
    set(colorbar_h,'YTickLabel',{sprintf('%5.2f',min_intensity_seg_all) sprintf('%5.2f',mean([min_intensity_seg_all max_intensity_seg_all])) ...
        sprintf('%5.2f',max_intensity_seg_all)});
    
end

h = [];
h_lines=[];

for i=1:length(chr_change)+1
    i
    if isempty(handles.subplots)
        h(i) = axes('Position',[left+0.02 bottom+0.01+(1/(length(chr_change)+2))*height*(i-1) width-left-0.02 (1/(length(chr_change)+2))*height*0.9]);
    else
        h(i) = handles.subplots(i);
        axes(h(i));
        cla reset
    end
    set(gca, 'Color', 'none');
    axis off
    if ~is_empty && isempty(handles.ChromsData{i})
        continue;
    end
    set(gca, 'ButtonDownFcn', 'mapAxesClick_Callback');
    %set(gca,'XTickLabel',[], 'YTickLabel',[]);
    %set(gca,'XTick',[], 'YTick',[]);
    set(gca, 'UserData', i);
%     axis off; 
    h_lines=[drawChrom(handles, i, 0, is_empty) h_lines];
end



