function EventsBinXL(EventData, currSlice, sliceCells, sliceFile)

    EventsBinHeaders = cell(1,numel(sliceCells));
    EventsBinValues  = cell(numel([EventData(sliceCells(1)).EventsPerBin.events]),numel(sliceCells));
    EventsBinTotals  = cell(1,numel(sliceCells));
            
    EventsBins = [{'Bins'};{EventData(sliceCells(1)).EventsPerBin.bins}'];
    TotalCells = 0;        
    for k = 1:numel(sliceCells)
        EventsBinHeaders{k}  = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).EventsPerBin)
            EventsBinValues(:,k) = {EventData(sliceCells(k)).EventsPerBin.events}';
            EventsBinTotals{k}  = sum([EventData(sliceCells(k)).EventsPerBin.events]);
            TotalCells = TotalCells + 1;
        else
            EventsBinValues(:,k) = num2cell(nan(numel({EventData(sliceCells(1)).EventsPerBin.bins}),1));
            EventsBinTotals{k}  = NaN;
        end
    end
            
    EventsBinOutput = horzcat(EventsBins, [EventsBinHeaders; EventsBinValues]);
    EventsBinTotals = [{'Total Events'}, EventsBinTotals];
    EventsBinOutput = [EventsBinOutput; EventsBinTotals];
            
    % Calculate means
    EventsBinMeans = nanmean(cell2mat(EventsBinValues),2);
    EventsBinTotalMean = nanmean(cell2mat(EventsBinTotals(2:end)));
    EventsBinMeansInclude = [{'Mean'}; num2cell(EventsBinMeans); {EventsBinTotalMean}];
    EventsBinOutput = [EventsBinOutput, EventsBinMeansInclude];
    
    % Calculate SE
    EventsBinSE = nanstd(cell2mat(EventsBinValues),0,2)./sqrt(TotalCells);
    EventsBinTotalSE = nanstd(cell2mat(EventsBinTotals(2:end)))./sqrt(TotalCells);
    EventsBinSEInclude = [{'SE'}; num2cell(EventsBinSE); {EventsBinTotalSE}];
    EventsBinOutput = [EventsBinOutput, EventsBinSEInclude];
    
    if ispc
        xlswrite(sliceFile, EventsBinOutput, 'EventsPerBin');
    else
        xlwrite(sliceFile, EventsBinOutput, 'EventsPerBin');
    end
end