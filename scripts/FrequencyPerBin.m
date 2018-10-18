function FrequencyPerBinOutput = FrequencyPerBin(currCell)

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
    % This is a counter for the for loop for indexing FrequencyPerBinOutput
    counter = 1;
    % Struct for containing Bin Names and FrequencyPerBin calculations
    FrequencyPerBinOutput = struct;

    % Cycle through for loop to determine Events Per Bin
    for i = startTime:binTime:endTime
        if i ~= endTime
            minTime = i;
            % Determine maximum sample index for current bin
            maxTime = i + binTime - samplingInterval;
            if maxTime > endTime; maxTime = endTime; end

            % Create bin name for current bin
            FrequencyPerBinOutput(counter).bins = strcat(num2str(minTime),'-',num2str(round(maxTime)));

            % Determine events in the bin (peak indices that fall between the
            % minIndex and maxIndex)
            EventsInBin = find(events >= minTime & events <= maxTime);


            % Calculate Frequencies
            % If there are no events in the bin, frequency is 0.
            if isempty(EventsInBin)
                FrequencyPerBinOutput(counter).frequency = NaN;
            % If there are events in the bin, calculate mean frequency.
            else
                % Extract periods of events in bin
                periodsInBin = diff(events(EventsInBin));
                % Calculate frequencies
                frequenciesInBin = 1./periodsInBin;
                % Calculate mean frequency
                if isempty(frequenciesInBin)
                    FrequencyPerBinOutput(counter).frequency = NaN;
                else
                    FrequencyPerBinOutput(counter).frequency = nanmean(frequenciesInBin);
                end
            end

            % Add 1 to counter
            counter = counter + 1;
        end
    end



end