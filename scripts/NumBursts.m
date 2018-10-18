function NumBurstsOutput = NumBursts(currCell)

burstIDs = unique(currCell.Bursts(:,1));
NumBurstsOutput = numel(burstIDs);


end