function h_lines=drawChrom(handles, chrom, standAlone, is_empty)

if ~is_empty
    Data=handles.ChromsData{chrom}.Data;
    if isfield(handles.ChromsData{chrom},'Genotypes')
        Genotypes=handles.ChromsData{chrom}.Genotypes;
    end
    if ~isempty(handles.PairedSample)
        pairedGenotypes=handles.PairedSample.DispStruct.Chrom{chrom}.Genotypes;
    end
    Locs=handles.ChromsData{chrom}.Locs;
    Segments=handles.ChromsData{chrom}.Segments;
    min_intensity_seg_all=handles.min_intensity_seg_all;
    max_intensity_seg_all=handles.max_intensity_seg_all;
end

%chr = handles.chr;
%location = handles.location;
% start_over = handles.start_over;
% end_over = handles.end_over;
% start_under = handles.start_under;
% end_under = handles.end_under;
chr_band_end = handles.chr_band_end;
% chr_band_start = handles.chr_band_start;
% chr_band_names = handles.chr_band_names;
chr_num = handles.chr_num;
end_p_location = handles.end_p_location;

ix=find(chr_num==chrom);
%chr_idx = find(chr==chrom);

% TOL = 0.000000000000001;
% qvals_over = handles.qvals_over;
% qvals_under = handles.qvals_under;
% qvals_over(qvals_over<TOL) = TOL;
% qvals_under(qvals_under<TOL) = TOL;
% plot_qvals_over = -log2(qvals_over);
% plot_qvals_under = -log2(qvals_under);
xmin = 0;
if is_empty || isempty(Data)
    ymin=0;
    ymax=1;
else
    ymin =  min(Data(:)); %min(min(plot_qvals_over(chr_idx)), min(plot_qvals_under(chr_idx)));
    if ymin>0
        ymin=0;
    end
    ymax =  max(Data(:)); %max(max(plot_qvals_over(chr_idx)), max(plot_qvals_under(chr_idx)));
    if ymax==0
        ymax=1;
    end
end

if standAlone
    xmax = max(chr_band_end(ix));
    set(handles.chromAxes, 'UserData', chrom);
    axes(handles.chromAxes);
    set(handles.chromSlider, 'Min', xmin, 'Max', xmax, 'Value', xmin);
else
    xmax = max(chr_band_end);
end

if  ymin==ymax
    ymax=ymin+0.0001;
end

axis_pos = [xmin xmax ymin ymax];

axis(axis_pos);
ax=axis;
hold on;
miny=ax(3);
maxy=ax(4);
lengthy = maxy-miny;

linewidth1=0.5;
if standAlone
    miny=miny-lengthy*0.125; %was 0.025
    lengthy=lengthy*1.25;    %was 1.05
    set(gca,'ylim',[miny miny+lengthy]); %was +-0.25%
    ax=axis; %for making band lines on entire chromosome height
    linewidth1=1.5;
end

if ~is_empty
    %over_idx = find(chr(start_over)==chrom);
    
    colormap('default');
    cmap = colormap;
    if(~isempty(Segments))
        if max_intensity_seg_all==min_intensity_seg_all
            color_segments=ones(size(Segments,1),1)*round(size(cmap,1)/2);
        else
            color_segments=ceil( (Segments(:,3)-min_intensity_seg_all) / (max_intensity_seg_all-min_intensity_seg_all) * (size(cmap,1)-1) )+1;
        end
        color_segments(color_segments>size(cmap,1))=size(cmap,1);        
    end
    
    for j=1:size(Segments,1)   %length(over_idx)
        rectangle_length = max(Segments(j,2) - Segments(j,1),1);  %location(end_over(over_idx(j))) - location(start_over(over_idx(j))) +1;
        rectangle('Position',[Segments(j,1),miny, rectangle_length,lengthy],'FaceColor',cmap(color_segments(j),:), 'EdgeColor', 'none' );
    end
    %under_idx = find(chr(start_under)==chrom);
%     for j=1:length(under_idx)
%         rectangle_length = location(end_under(under_idx(j))) - location(start_under(under_idx(j))) +1;
%         rectangle('Position',[location(start_under(under_idx(j))),miny, rectangle_length,lengthy],'FaceColor',[0.51 0.6 0.98], 'EdgeColor', 'none' );
%     end
    
    hold on;
    if strcmp(handles.SampleName,'All Samples')
        lines_colormap=[0 0 0; 1 0 0];%black - num of samples involved in deletion. red - amplification.
    else
        lines_colormap=colormap('lines');
    end
    
    Locs_idx_p=find(Locs <= end_p_location(chrom));
    Locs_idx_q=find(Locs > end_p_location(chrom));
    
    lines_gap_p=ones(length(Locs_idx_p),1)*lengthy/50;
    lines_gap_q=ones(length(Locs_idx_q),1)*lengthy/50;
    
    for j=1:size(Data,2)
        if(j==2) % only the two chromosomes copy numbers should be separated
            plot(Locs(Locs_idx_p), Data(Locs_idx_p,j)+lines_gap_p*(j-1), 'color',lines_colormap(j,:),'linewidth',linewidth1);
            plot(Locs(Locs_idx_q), Data(Locs_idx_q,j)+lines_gap_q*(j-1), 'color',lines_colormap(j,:),'linewidth',linewidth1);
        else
            plot(Locs(Locs_idx_p), Data(Locs_idx_p,j), 'color',lines_colormap(j,:),'linewidth',linewidth1);
            plot(Locs(Locs_idx_q), Data(Locs_idx_q,j), 'color',lines_colormap(j,:),'linewidth',linewidth1);
        end
    end
    %plot(location(chr_idx), plot_qvals_under(chr_idx), 'k');
    colormap('default');
end

rectangle('Position',[1,miny, end_p_location(chrom),lengthy],'Curvature',0.4, 'LineWidth',2,'LineStyle','-');
rectangle('Position',[end_p_location(chrom)+1,miny, max(chr_band_end(ix))-end_p_location(chrom),lengthy],'Curvature',0.4, 'LineWidth',2,'LineStyle','-');

h_lines=[];
for j=1:length(ix)-1
    if isequal(chr_band_end(ix(j)), end_p_location(chrom))
        continue;
    end
    h_lines(end+1)=line([chr_band_end(ix(j)) chr_band_end(ix(j))], [ax(3) ax(4)],'LineWidth',1,'Color', 'k');
end

if chrom<23
    chrLabel = num2str(chrom);
elseif chrom==23
    chrLabel = 'X';
else
    chrLabel = 'Y';
end

if standAlone
    %plot genotype
    if exist('Genotypes','var')
        idx=find(Genotypes==0 | Genotypes==3);
        AA_BB_perc=num2str(round(length(idx)/length(Genotypes)*100));
        %plot(Locs(idx), lengthy*0.9, '*','color','g');
        scatter(Locs(idx), ones(1,length(idx))*lengthy*0.89,25,'g','.'); %AA/BB
        idx=find(Genotypes==1 | Genotypes==2);
        AB_perc=num2str(round(length(idx)/length(Genotypes)*100));
        scatter(Locs(idx), ones(1,length(idx))*lengthy*0.88,25,'r','.'); %AB
        
        if ~isempty(handles.PairedSample)
            idx=find((Genotypes==0 | Genotypes==3) & (pairedGenotypes==1 | pairedGenotypes==2));
            LOH_perc=num2str(round(length(idx)/length(Genotypes)*100));
            scatter(Locs(idx), ones(1,length(idx))*lengthy*0.9,25,'m','.'); %LOH
        end
    end
    %plot(Locs(idx), lengthy*0.9, '*','color','r'); %plot causes bug.
    
%     ylim1=get(gca,'ylim');
%     delta_y=(ylim1(2) - ylim1(1)) * 0.05;
%     ylim1=ylim1+[-delta_y delta_y];
%     set(gca,'ylim',ylim1);
%    axis tight
    if strcmp(handles.SampleName,'All Samples')
        ylabel('Num Del/Amp Samples');
    else
        ylabel('-log2(q-value)');
    end
    title1=['Chromosome   ' chrLabel];
    if exist('Genotypes','var')
        title1=[title1 '       \color{green}*AA/BB(' AA_BB_perc '%)      \color{red}*AB(' AB_perc '%)'];
    end
    if ~isempty(handles.PairedSample)
        title1=[title1 '       \color{magenta}*LOH(' LOH_perc '%)'];
    end
    title(title1);
    callbackFunction = 'chromAxesClick_Callback';
    set(handles.chromSlider, 'Min', 0, 'Max', xmax);
    updateChromView(handles);
else
    text(-5000000, miny+(lengthy/2), chrLabel, 'FontSize', 12, 'FontWeight', 'bold');
    callbackFunction = 'mapAxesClick_Callback';
end

h=get(gca, 'Children');
set(h, 'ButtonDownFcn', callbackFunction);

% if standAlone
%     cmenu = uicontextmenu;
%     uimenu(cmenu, 'Label', 'Export', 'Callback', 'export_Callback');
%     set(h,'UIContextMenu', cmenu);
%     set(gca,'UIContextMenu', cmenu);
% end

