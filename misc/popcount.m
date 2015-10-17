% Count number of bits in each element of X 
function PC = popcount(X)

PC = 0;
for i=1:32
    PC = PC + bitget(X,i);
end
