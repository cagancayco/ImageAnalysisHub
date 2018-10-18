function InterBurstIntervalOutput = InterBurstInterval(currCell)
    [IDs, idx, ~] = unique(currCell.Bursts(:,1));
    InterBurstIntervalOutput = zeros(numel(IDs)-1,1);
    
    for i = 1:numel(IDs)-1
        burstInterPeriod = currCell.Bursts(idx(i+1),2) - currCell.Bursts(idx(i+1)-1,2);
        InterBurstIntervalOutput(i) = burstInterPeriod;
    end


end