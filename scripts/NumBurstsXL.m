function NumBurstsXL(EventData, currSlice, sliceCells, sliceFile)

    NumBurstsHeaders = cell(1,numel(sliceCells));
    NumBurstsValues = cell(1,numel(sliceCells));
    
    TotalCells = 0;
    for k = 1:numel(sliceCells)
        NumBurstsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).nBursts)
            NumBurstsValues{1,k} = EventData(sliceCells(k)).nBursts;
            TotalCells = TotalCells + 1;
        else
            NumBurstsValues{1,k} = NaN;
        end
    end

    NumBurstsOutput = vertcat(NumBurstsHeaders,NumBurstsValues);
    
    % Calculate mean
    NumBurstsMean = nanmean(cell2mat(NumBurstsValues));
    NumBurstsOutput = [NumBurstsOutput,[{'Mean'};num2cell(NumBurstsMean)]];
    
    % Calculate SE
    NumBurstsSE = nanstd(cell2mat(NumBurstsValues))./sqrt(TotalCells);
    NumBurstsOutput = [NumBurstsOutput, [{'SE'};num2cell(NumBurstsSE)]];
    
    if ispc
        xlswrite(sliceFile, NumBurstsOutput, 'NumBursts');
    else
        xlwrite(sliceFile, NumBurstsOutput, 'NumBursts');
    end
    
    
end