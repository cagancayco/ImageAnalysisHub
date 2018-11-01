function FractionSingleEventsXL(EventData, currSlice, sliceCells, sliceFile)

    FracSingEventsHeaders = cell(1,numel(sliceCells));
    FracSingEventsValues  = cell(1,numel(sliceCells));
    
    TotalCells = 0;
    for k = 1:numel(sliceCells)
        FracSingEventsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).fracSingleEvents)
            FracSingEventsValues{1,k} = EventData(sliceCells(k)).fracSingleEvents;
            TotalCells = TotalCells + 1;
        else
            FracSingEventsValues{1,k} = NaN;
        end
    end

    FracSingEventsOutput = vertcat(FracSingEventsHeaders,FracSingEventsValues);
    
    % Calculate mean
    FracSingEventsMean = nanmean(cell2mat(FracSingEventsValues));
    FracSingEventsOutput = [FracSingEventsOutput, [{'Mean'};num2cell(FracSingEventsMean)]];
    
    % Calculate SE
    FracSingEventsSE = nanstd(cell2mat(FracSingEventsValues))./sqrt(TotalCells);
    FracSingEventsOutput = [FracSingEventsOutput, [{'SE'}; num2cell(FracSingEventsSE)]];
    
    FracSingEventsValues = FracSingEventsValues';
    emptyCol = cell(size(FracSingEventsValues));
    emptyRows = cell(size(FracSingEventsValues,1)-size(FracSingEventsOutput,1),size(FracSingEventsOutput,2));
    
    FracSingEventsOutput = [[FracSingEventsOutput;emptyRows],emptyCol,FracSingEventsValues];
    
    
    if ispc
        xlswrite(sliceFile, FracSingEventsOutput, 'FracSingleEvents');
    else
        xlwrite(sliceFile, FracSingEventsOutput, 'FracSingleEvents');
    end
end