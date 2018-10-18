function BurstFrequencyOutput = BurstFrequency(currCell)

    BurstFrequencyOutput = struct;
    
    BurstFrequencyOutput.inter = struct;
    [IDs, idx, ~] = unique(currCell.Bursts(:,1));

    BurstFrequencyOutput.intra.burstID = IDs;
    
    BurstFrequencyOutput.intra.frequency = zeros(numel(IDs),1);
    %BurstFrequencyOutput.inter = zeros(numel(IDs)-1,1);
    for i = 1:numel(IDs)
        if i ~= numel(IDs)
            burstEvents = currCell.Bursts(idx(i):(idx(i+1)-1),2);
            %burstInterPeriod = currCell.Bursts(idx(i+1),2) - currCell.Bursts(idx(i+1)-1,2);
            %BurstFrequencyOutput.inter(i) = 1/burstInterPeriod;
            burstDur = currCell.Bursts(idx(i+1)-1,2) - currCell.Bursts(idx(i),2);
        else
            burstEvents = currCell.Bursts(idx(i):end,2);
            burstDur = currCell.Bursts(end,2) - currCell.Bursts(idx(i),2);
        end
        %burstIntraPeriods = diff(burstEvents);
        %burstIntraFrequencies = 1./burstIntraPeriods;
        nEvents = numel(burstEvents);
        BurstFrequencyOutput.intra.frequency(i) = nEvents/burstDur;
    end

    


end