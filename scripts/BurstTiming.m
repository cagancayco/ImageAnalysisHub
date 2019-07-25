function BurstTimingOutput = BurstTiming(currCell)

[~, idx, ~] = unique(currCell.Bursts(:,1));

BurstTimingOutput = currCell.Bursts(idx, 2);


end