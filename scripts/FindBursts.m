function BurstsOutput = FindBursts(currCell)

    if strcmp(currCell.algorithm,'minEASE')
        events = currCell.events;
        events = currCell.time(events)';
        
        classesToInc = currCell.classesToInc;
        eventClass   = currCell.eventClass;
        eventsToInc = zeros(length(events),1);
        for i = 1:length(eventsToInc)
            if ~isempty(find(classesToInc == eventClass(i),1))
                eventsToInc(i) = 1;
            end
        end
        eventsToInc = logical(eventsToInc);
        
        events = events(eventsToInc);
    else
        events = currCell.events;
    end

    startTime = currCell.startTime;
    endTime   = currCell.time(end);

    if currCell.endTime ~= 0   % If endTime is 0; then endTimeIndex is last index in data and timepoints vectors
        endTime = currCell.endTime; % end time in seconds
    end
    
    events(events < startTime) = [];
    events(events > endTime)   = [];
    
    if numel(events) >= 3
        BurstsOutput = zeros(numel(events),2);
        currID = 1;

        BurstsOutput(1,1) = currID;
        BurstsOutput(1,2) = events(1);

        for i = 2:length(BurstsOutput)

            if abs(events(i) - events(i-1)) <= currCell.threshold
                BurstsOutput(i,1) = currID;
                BurstsOutput(i,2) = events(i);
            else
                currID = currID + 1;
                BurstsOutput(i,1) = currID;
                BurstsOutput(i,2) = events(i);
            end

        end

        IDs = unique(BurstsOutput(:,1));
        IDcounts = histc(BurstsOutput(:,1),1:numel(IDs));

        rmBurst = find(IDcounts < currCell.burstThreshold);
        rmEvents = ismember(BurstsOutput(:,1),IDs(rmBurst));
        BurstsOutput(rmEvents,:) = [];

        [~,~,newIDs] = unique(BurstsOutput(:,1));
        BurstsOutput(:,1) = newIDs;

    else
        BurstsOutput = [];
    end

end