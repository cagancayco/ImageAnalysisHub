function BurstDurationOutput = BurstDuration(currCell)

    [IDs,idx,~] = unique(currCell.Bursts(:,1));
    
    BurstDurationOutput = struct;
    BurstDurationOutput.burstID = IDs;
    
    BurstDurationOutput.duration = zeros(numel(IDs),1);
    for i = 1:numel(IDs)
        if i ~= numel(IDs)
            BurstDurationOutput.duration(i) = currCell.Bursts(idx(i+1)-1,2) - currCell.Bursts(idx(i),2);
        else
            BurstDurationOutput.duration(i) = currCell.Bursts(end,2) - currCell.Bursts(idx(i),2);
        end
    end


end