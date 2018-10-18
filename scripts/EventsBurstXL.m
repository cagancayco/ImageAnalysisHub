function EventsBurstXL(EventData, currSlice, sliceCells, sliceFile)

    EventsBurstHeaders = cell(1,numel(sliceCells));
    EventsBurstResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    EventsBurstStats = cell(5, numel(sliceCells) + 2);
    EventsBurstStats{1,1} = 'mean';
    EventsBurstStats{2,1} = 'SE';
    EventsBurstStats{3,1} = 'median';
    EventsBurstStats{4,1} = 'max';
    EventsBurstStats{5,1} = 'min';
    
    for k = 1:numel(sliceCells)
        EventsBurstHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).nEventsPerBurst)
            EventsBurstResults{k} = [EventData(sliceCells(k)).nEventsPerBurst.nEvents];
            rowCounts(k) = numel([EventData(sliceCells(k)).nEventsPerBurst.nEvents]);
            
            % Mean
            EventsBurstStats{1,k+1} = mean(EventsBurstResults{k});
            % SE
            EventsBurstStats{2,k+1} = std(EventsBurstResults{k})./sqrt(rowCounts(k));
            % Median
            EventsBurstStats{3,k+1} = median(EventsBurstResults{k});
            % Max
            EventsBurstStats{4,k+1} = max(EventsBurstResults{k});
            % Min
            EventsBurstStats{5,k+1} = min(EventsBurstResults{k});
        else
            rowCounts(k) = 0;
        end 
    end

    EventsBurstOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        EventsBurstOutput(1:rowCounts(k),k) = num2cell(EventsBurstResults{k});
    end
    TotalEventsBurst = reshape(EventsBurstOutput, [numel(EventsBurstOutput),1]);
    TotalEventsBurst = reshape(cell2mat(TotalEventsBurst), [sum(rowCounts),1]);
    EventsBurstOutput(:,end) = num2cell(TotalEventsBurst);
    
    % Mean
    EventsBurstStats{1,end} = mean(TotalEventsBurst);
    % SE
    EventsBurstStats{2,end} = std(TotalEventsBurst)./sqrt(numel(TotalEventsBurst));
    % Median
    EventsBurstStats{3,end} = median(TotalEventsBurst);
    % Max
    EventsBurstStats{4,end} = max(TotalEventsBurst);
    % Min
    EventsBurstStats{5,end} = min(TotalEventsBurst);
    
    EventsBurstOutput = [[{'Descriptive Stats'}, EventsBurstHeaders, {'All Events Per Burst'}]; EventsBurstStats; [{''}, EventsBurstHeaders, {'All Events Per Burst'}]; [cell(numel(TotalEventsBurst),1), EventsBurstOutput]];

    if ispc
        xlswrite(sliceFile, EventsBurstOutput, 'EventsPerBurst');
    else
        xlwrite(sliceFile, EventsBurstOutput, 'EventsPerBurst');
    end
end