% Read a binary file 
function dat_s = ReadBin(file_s)
fid=fopen(file_s,'rb');
dat_s=fread(fid,'float32');
fclose(fid);
