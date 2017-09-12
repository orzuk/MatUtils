% Determine width and length we should take when showing n figs
function [w l] = subplot_dims(n)
w = ceil(sqrt(n));
if(w*(w-1) < n)
    l = w;
else
    l = w-1;
end
