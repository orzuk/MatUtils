% Unite several GWAS studies of the same disease. 
% Remove duplicate loci (another alternative: keep only largest study!)
% 
% Input: 
% data - structure with data for each SNP
% 
% Output: 
% data_united - unite all SNPs for the same disease
%
function data_united = unite_gwas_data_by_traits(data) 

[unique_traits unique_inds] = unique(data.Trait); % Find manually which ones are quantitative or qualitative traits

% data_united = unique_struct_by_field(data, 'Disease_Trait'); 

bad_inds = []; 
for i=1:length(unique_inds)
	cur_trait_inds = strmatch(unique_traits{i}, data.Trait,  'exact');     
	cur_regions = data.Region(cur_trait_inds);
%	cur_pubmedid = data.PUBMEDID(cur_trait_inds); 
	
	% Pick one SNP from each region (too complicated to do otherwise). Pick the one with highest OR * MAF (?)
	cur_snp_power = data.OR(cur_trait_inds) .* ...
        min(data.MAF(cur_trait_inds), 1-data.MAF(cur_trait_inds));
		
	unique_regions = unique(cur_regions); 
	for j=1:length(unique_regions)
		cur_region_inds = strmatch(unique_regions{j}, cur_regions, 'exact');
		[~, max_snp_in_region] = max(cur_snp_power(cur_region_inds));			
		cur_bad_inds = cur_region_inds(setdiff(1:length(cur_region_inds), max_snp_in_region)); 
		bad_inds = [bad_inds vec2row(cur_trait_inds(cur_bad_inds))];
	end
end
good_inds = setdiff(1:length(data.Trait), bad_inds); 
data_united = struct_by_inds(data, good_inds); 

got_totally = length(good_inds)
