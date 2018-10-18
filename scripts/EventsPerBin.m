function EventsPerBinOutput = EventsPerBin(currCell)

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
    

    startTime = 0;
    endTime   = currCell.time(end);

    if currCell.endTime ~= 0   % If endTime is 0; then endTimeIndex is last index in data and timepoints vectors
        endTime = currCell.endTime; % end time in seconds
    end  
    
    binTime = currCell.binTime;
    samplingInterval = currCell.sampInt;
    % This is a counter for the for loop for indexing EventsPerBinOutput
    counter = 1;
    % Struct for containing Bin Names and EventsPerBin calculations
    EventsPerBinOutput = struct;


    % Cycle through for loop to determine Events Per Bin
    for i = startTime:binTime:endTime
        if i ~= endTime
            minTime = i;
            % Determine maximum sample index for current bin
            maxTime = i + binTime - samplingInterval;
            if maxTime > endTime; maxTime = endTime; end

            % Create bin name for current bin
            EventsPerBinOutput(counter).bins = strcat(num2str(minTime),'-',num2str(round(maxTime)));

            % Calculate number of events in the bin (number of peak indices that
            % fall between minIndex and maxIndex)
            numEventsInBin = numel(events(events >= minTime & events <= maxTime));

            % Fill in output vector
            EventsPerBinOutput(counter).events = numEventsInBin;

            % Add 1 to counter
            counter = counter + 1;
        end
    end


end