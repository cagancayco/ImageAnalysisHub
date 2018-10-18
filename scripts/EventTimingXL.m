function EventTimingXL(EventData, currSlice, sliceCells, sliceFile)

    EventTimingHeaders = cell(1,numel(sliceCells));
    EventTiming = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
            
    for k = 1:numel(sliceCells)
        EventTimingHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).EventTiming)
            EventTiming{k} = [EventData(sliceCells(k)).EventTiming];
            rowCounts(k)  = numel([EventData(sliceCells(k)).EventTiming]);
        else
            rowCounts(k) = 0;
        end
    end
            
    EventTimingOutput = cell(max(rowCounts),numel(sliceCells));
    for k = 1:numel(sliceCells)
        EventTimingOutput(1:rowCounts(k),k) = num2cell(EventTiming{k});
    end
            
    EventTimingOutput = [EventTimingHeaders; EventTimingOutput];
    
    if ispc
        xlswrite(sliceFile, EventTimingOutput, 'Event Timing');
    else
        xlwrite(sliceFile, EventTimingOutput, 'Event Timing');
    end        
end