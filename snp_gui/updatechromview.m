function updateChromView(handles)

chrom = get(handles.chromAxes,'UserData');
if isempty(chrom)
    switch get(get(handles.uipanel_choose_view,'SelectedObject'), 'Tag')   % Get Tag of selected object
        case 'radioNone'
            lines_on_off(handles,'off');
        otherwise
            lines_on_off(handles,'on');
    end
    return;
end
    

ix=find(handles.chr_num==chrom);

%xmax = max(handles.chr_band_end(ix));
ax=axis;
lastSelectedTag = get(handles.uipanel_choose_view,'UserData');
switch get(get(handles.uipanel_choose_view,'SelectedObject'), 'Tag')   % Get Tag of selected object
    case 'radioGenes'
        delete(findobj(gcf,'Tag',lastSelectedTag)); 
        chr_idx = find(handles.annot.chr==chrom);
        current_idx = find(handles.annot.start(chr_idx)<=ax(2) & handles.annot.end(chr_idx)>=ax(1));
        x_middle = handles.annot.start(chr_idx(current_idx)) + (handles.annot.end(chr_idx(current_idx)) - handles.annot.start(chr_idx(current_idx)))/2;
        [sorted,ix]=unique(x_middle);
        set(gca,'xtick',sorted,'XTickLabel','');
        text(sorted, repmat(ax(3)-0.17*(ax(4)-ax(3)), length(ix),1), handles.annot.symbols(chr_idx(current_idx(ix))), 'Rotation', -90, 'FontSize', 7, 'Tag', 'geneLabels');
        lastSelectedTag = 'geneLabels';
        lines_on_off(handles,'on');
    case 'radioBands'
        delete(findobj(gcf,'Tag',lastSelectedTag)); 
        set(gca,'xtick',[],'XTickLabel','');
        text(handles.chr_band_start(ix) + (handles.chr_band_end(ix) - handles.chr_band_start(ix))/2, repmat(ax(3)-0.17*(ax(4)-ax(3)),length(ix),1), handles.chr_band_names(ix), 'Rotation', -90, 'FontSize', 8, 'Tag', 'bandLabels');
        lastSelectedTag = 'bandLabels';
        lines_on_off(handles,'on');
    case 'radioNone'
        delete(findobj(gcf,'Tag',lastSelectedTag)); 
        set(handles.chromAxes,'XTickMode','auto');
        xTickNum = get(handles.chromAxes, 'xTick');
        text(xTickNum', repmat(ax(3)-0.17*(ax(4)-ax(3)),length(xTickNum),1), num2str(xTickNum'), 'Rotation', -90, 'FontSize', 8, 'Tag', 'Location');
        lastSelectedTag = 'Location';
        lines_on_off(handles,'off');
    case 'radioSNP'
        delete(findobj(gcf,'Tag',lastSelectedTag)); 
        chr_idx = find(handles.chip_annot.chr_vec==chrom);
        current_idx = find(handles.chip_annot.chr_loc_vec(chr_idx)<=ax(2) & handles.chip_annot.chr_loc_vec(chr_idx)>=ax(1));
        set(gca,'xtick',handles.chip_annot.chr_loc_vec(chr_idx(current_idx)),'XTickLabel','');
        text(handles.chip_annot.chr_loc_vec(chr_idx(current_idx)), repmat(ax(3)-0.17*(ax(4)-ax(3)),length(chr_idx(current_idx)),1), handles.chip_annot.snps_gene_symbols(chr_idx(current_idx)), 'Rotation', -90, 'FontSize', 7, 'Tag', 'SNP');
        lastSelectedTag = 'SNP';
        lines_on_off(handles,'on');
end
set(handles.uipanel_choose_view, 'UserData', lastSelectedTag);

guidata(gcbo, handles);

function lines_on_off(handles,on_or_off)

if ~strcmp(get(handles.h_lines(1),'visible'),on_or_off)
    for i=1:length(handles.h_lines)
        set(handles.h_lines(i),'visible',on_or_off);
    end
end
    
if isfield(handles,'h_lines_standAlone') && ishandle(handles.h_lines_standAlone(1)) && ...
        ~strcmp(get(handles.h_lines_standAlone(1),'visible'),on_or_off)
    for i=1:length(handles.h_lines_standAlone)
        set(handles.h_lines_standAlone(i),'visible',on_or_off);
    end
end


