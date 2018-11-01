function TimeFirstBurstXL(EventData, currSlice, sliceCells, sliceFile)

    TimeFirstBurstHeaders = cell(1, numel(sliceCells));
    TimeFirstBurstValues = cell(1, numel(sliceCells));
    
    TotalCells = 0;
    for k = 1:numel(sliceCells)
        TimeFirstBurstHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).timeFirstBurst)
            TimeFirstBurstValues{1,k} = EventData(sliceCells(k)).timeFirstBurst;
            TotalCells = TotalCells + 1;
        else
            TimeFirstBurstValues{1,k} = NaN;
        end
    end

    TimeFirstBurstOutput = vertcat(TimeFirstBurstHeaders,TimeFirstBurstValues);
    
    % Calculate mean
    TimeFirstBurstMean = nanmean(cell2mat(TimeFirstBurstValues));
    TimeFirstBurstOutput = [TimeFirstBurstOutput,[{'Mean'};num2cell(TimeFirstBurstMean)]];
    
    % Calculate SE
    TimeFirstBurstSE = nanstd(cell2mat(TimeFirstBurstValues))./sqrt(TotalCells);
    TimeFirstBurstOutput = [TimeFirstBurstOutput, [{'SE'};num2cell(TimeFirstBurstSE)]];
  
    TimeFirstBurstValues = TimeFirstBurstValues';
    emptyCol = cell(size(TimeFirstBurstValues));
    emptyRows = cell(size(TimeFirstBurstValues,1)-size(TimeFirstBurstOutput,1),size(TimeFirstBurstOutput,2));
    
    TimeFirstBurstOutput = [[TimeFirstBurstOutput;emptyRows],emptyCol,TimeFirstBurstValues];
    
    if ispc
        xlswrite(sliceFile, TimeFirstBurstOutput, 'TimeToFirstBurst');
    else
        xlwrite(sliceFile, TimeFirstBurstOutput, 'TimeToFirstBurst');
    end
end