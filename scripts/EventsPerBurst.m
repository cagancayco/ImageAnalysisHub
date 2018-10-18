function EventsPerBurstOutput = EventsPerBurst(currCell)

    IDs = unique(currCell.Bursts(:,1));
    
    EventsPerBurstOutput = struct;
    EventsPerBurstOutput.burstID = IDs;
    EventsPerBurstOutput.nEvents = histc(currCell.Bursts(:,1), 1:numel(IDs));



end