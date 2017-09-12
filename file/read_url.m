% Read a file given as an html url
function lines = read_url(url_str)

url = java.net.URL(url_str);

is = openStream(url);
isr = java.io.InputStreamReader(is);
br = java.io.BufferedReader(isr);

lines = {}; i=1; lines{1} = char(readLine(br));
while(~isempty(lines{i}))
    i=i+1; lines{i} = char(readLine(br));
end

