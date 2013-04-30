function hapmap_data_http = hapmap_version_to_http(hapmap_version)

% Choose correct url
if(strcmp(hapmap_version, 'r21'))
    hapmap_data_http = 'http://www.hapmap.org/genotypes/2006-07/fwd_strand/non-redundant/';
end
if(strcmp(hapmap_version, 'r21a'))
    hapmap_data_http = 'http://www.hapmap.org/genotypes/latest_ncbi_build35/fwd_strand/non-redundant/';
end
if(strcmp(hapmap_version, 'r22'))
    hapmap_data_http = 'http://www.hapmap.org/genotypes/latest_ncbi_build36/fwd_strand/non-redundant/';
end
