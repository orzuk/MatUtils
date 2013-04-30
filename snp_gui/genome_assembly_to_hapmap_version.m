% return the current hapmap version 
function hapmap_version = genome_assembly_to_hapmap_version()

hapmap_version = 'ERR_unknown_ver';

genome_assembly = get_genome_assembly();
if(strmatch(genome_assembly, 'hg17'))
    hapmap_version = 'r21';
end    
if(strmatch(genome_assembly, 'hg18'))
    hapmap_version = 'r22';
end    

