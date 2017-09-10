% Get the file names of hapmap genotypes from the web
function file_names = get_hapmap_file_names_from_web(hapmap_version)

% Choose correct url
hapmap_data_http = hapmap_version_to_http(hapmap_version);

s = urlread(hapmap_data_http); % Read data from project
start_inds = strfind(s,'genotypes_chr'); end_inds = strfind(s,'txt.gz')+5;
if(size(start_inds) ~= size(end_inds))
    sprintf('Error In Page - Go and Download it Manually ...')
    file_names = [];
    return;
end

file_names = {};
for i=1:length(start_inds)
    file_names{i} = s(start_inds(i):end_inds(i));
end
file_names = unique(file_names);

if(isempty(strfind(file_names{1},hapmap_version)))
    sprintf('Error In Hapmap Version - Go and Download it Manually ...')
    file_names = [];
    return;
end
