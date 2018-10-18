function FrequencyBinXL(EventData, currSlice, sliceCells, sliceFile)
    FreqBinHeaders = cell(1,numel(sliceCells));
    FreqBinValues  = cell(numel([EventData(sliceCells(1)).FrequencyPerBin.frequency]),numel(sliceCells));
    FreqBinTotals  = cell(1,numel(sliceCells));
            
    FreqBins = [{'Bins'};{EventData(sliceCells(1)).FrequencyPerBin.bins}'];
    TotalCells = 0; 
    for k = 1:numel(sliceCells)
        FreqBinHeaders{k}  = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).FrequencyPerBin)
            FreqBinValues(:,k) = {EventData(sliceCells(k)).FrequencyPerBin.frequency}';
            FreqBinTotals{k}  = mean([EventData(sliceCells(k)).FrequencyPerBin.frequency]);
            TotalCells = TotalCells + 1;
        else
            FreqBinValues(:,k) = num2cell(nan(numel({EventData(sliceCells(1)).FrequencyPerBin.bins}),1));
            FreqBinTotals{k}  = NaN;
        end
    end
            
    FreqBinOutput = horzcat(FreqBins, [FreqBinHeaders; FreqBinValues]);
    FreqBinTotals = [{'Mean Frequency Per Cell'}, FreqBinTotals];
    FreqBinOutput = [FreqBinOutput; FreqBinTotals];
            
    % Calculate means
    FreqBinMeans = nanmean(cell2mat(FreqBinValues),2);
    FreqBinTotalMean = nanmean(cell2mat(FreqBinTotals(2:end)));
    FreqBinMeansInclude = [{'Mean Frequency Including Cells with no Events'}; num2cell(FreqBinMeans); {FreqBinTotalMean}];
    FreqBinOutput = [FreqBinOutput, FreqBinMeansInclude];
            
    % Calculate SE
    FreqBinSE = nanstd(cell2mat(FreqBinValues),0,2)./sqrt(TotalCells);
    FreqBinTotalSE = nanstd(cell2mat(FreqBinTotals(2:end)))./sqrt(TotalCells);
    FreqBinSEInclude = [{'SE'}; num2cell(FreqBinSE); {FreqBinTotalSE}];
    FreqBinOutput = [FreqBinOutput, FreqBinSEInclude];
    
    if ispc
        xlswrite(sliceFile, FreqBinOutput, 'FrequencyPerBin'); 
    else
        xlwrite(sliceFile, FreqBinOutput, 'FrequencyPerBin'); 
    end

end