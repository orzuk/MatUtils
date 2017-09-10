% Write a link in html format
%
% Input: 
% link_name - text appearing
% link_url - url linked
% 
% Output: 
% link_str - string containing link
% 
function link_str = html_write_link(link_name, link_url)

link_str = ['<a href="' link_url '">' link_name '</a>']; 
