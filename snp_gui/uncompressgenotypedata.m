% Move from bit representation to full mat
function UncompressedSnpsMat = UncompressGenotypeData(SnpsData)
UncompressedSnpsMat = zeros(size(SnpsData{1},1),180);
for i=1:3
    for j=1:30
        UncompressedSnpsMat(:,(i-1)*30+j) = bitget(SnpsData(:,i), j);
        UncompressedSnpsMat(:,90+(i-1)*30+j) = bitget(SnpsData(:,i+3), j);
        
    end
end

