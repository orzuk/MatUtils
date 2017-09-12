% Add an image in html format
function img_link_str = html_write_img(img_name, img_link, img_width, img_height)

if(~exist('img_width', 'var') || isempty(img_width))
    img_width = 500;
end
if(~exist('img_height', 'var') || isempty(img_height))
    img_height = 400;
end

img_link_str = [img_name ' <br> ' ...
    '<a href="' img_link '" target="_blank">' ...
    '<img src="' img_link '" ' ...
    'alt="' img_name '"width="' num2str(img_width) '" height="' num2str(img_height) '"></img> </a><br>']; % Add figure showing effect size explained

