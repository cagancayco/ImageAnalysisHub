function EventTimesOutput = EventTimes(currCell)

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


    % Find EventTimes by indexing timepoints vector with peakIndex
    EventTimesOutput = events(events >= startTime & events <= endTime);

end