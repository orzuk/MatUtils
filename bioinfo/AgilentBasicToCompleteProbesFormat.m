% Converts a basic probes format file to a complete format for agilent eArray. 
% The complete format gives you more information. 
% We can also filter probes based on quality scores etc.
%
% Input: 
% basic_probes_file - a file with probes in a basic format: probe-is sequence
%
% The output: 
% Written to complete_probes_outfile - a file with complete format: 
% ProbeID Sequence TargetID Accessions GeneSymbols Description ChromosomalLocation
%
function AgilentBasicToCompleteProbesFormat(basic_probes_file, complete_probes_outfile, genome_version)

% R = textread(basic_probes_file, '%s', 'delimiter', '\n'); % load input file (tab-delimited)
R = loadcellfile(basic_probes_file);
R = R(2:end,:); 
num_probes = size(R, 1); % get total number of probes 

T = cell(num_probes, 7); % set the cell file to save 

for i=1:num_probes
%     if(mod(i, 500) == 0)
%         i_is = i
%     end
    w = R{i,1};
    strand_ind = strfind(w, 'strand'); chr_ind = strfind(w, 'chr'); dash_ind = strfind(w, '-'); colon_ind = strfind(w, ':');
    if(isempty(dash_ind) || (max(dash_ind) < chr_ind))
        dash_ind = strfind(w, '.'); dash_ind = min(dash_ind(dash_ind > chr_ind)); % take the first dot 
    end
    T{i,1} = [R{i,1}(1:strand_ind-1) num2str(mod_max(i, 4)) '_' num2str(i) '_' genome_version]; % get the ProbeID (include genome_version to be unique)
    T{i,2} = R{i,3}; % get the sequence
    T{i,3} = R{i,1}(1:strand_ind-2); % target id
    T{i,4} = 'NA|00000'; %%% R{i,1}(1:strand_ind-2);  % Accessions (?) - we don't have one so just repeat target id
    T{i,5} = R{i,1}(1:strand_ind-2); % Gene symbols (?) - we don't have anything here so just put again target id 
    T{i,6} = [w(strand_ind:chr_ind-1) R{i,15} '_bp-start_' num2str(R{i,2}) '_bp-end_' num2str(R{i,5})]; % description (include strand , bc score and possible cross-hybs)
    if(~isempty(R{i,14}))
        T{i,6} = [T{i,6} '_X-Hyb_' R{i,14}];
    end
    T{i,6} = ['"' T{i,6} '"']; % add " 
    if(~isempty(colon_ind))
        T{i,7} = [w(chr_ind:dash_ind(end)-1) '-' num2str( str2num(w(colon_ind+1:dash_ind(end)-1)) + R{i,2} - 1)]; % chromosomal location
    else % whenever we cannot obtain chromosomal location
        T{i,7} = 'chr1:1-1';
    end
end


savecellfile(T, complete_probes_outfile); % save complete format (tab-delimited)


