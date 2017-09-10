% Returns 0 if false, 1 if true
function isposdef = IsPositiveDefinite(M)
isposdef = true;
for i=1:length(M)
    if ( det( M(1:i, 1:i) ) <= 0 )
        isposdef = false;
        break;
    end
end
