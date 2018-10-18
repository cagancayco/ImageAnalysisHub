function FractionalATFirstLastXL(EventData, currSlice, sliceCells, sliceFile)

    FracATFirstLastHeaders = cell(1,numel(sliceCells)+2);
    FracATFirstLastOutput  = cell(9,numel(sliceCells)+2);
    TotalCells = 0;
    
    for k = 1:numel(sliceCells)
        FracATFirstLastHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).FracATFirstLast)
            FracATFirstLastOutput{1,k} = EventData(sliceCells(k)).FracATFirstLast.nQuiescentPeriods;
            FracATFirstLastOutput{3,k} = EventData(sliceCells(k)).FracATFirstLast.sumQuiescentTime;
            FracATFirstLastOutput{5,k} = EventData(sliceCells(k)).FracATFirstLast.fracQuiescentTime;
            FracATFirstLastOutput{7,k} = EventData(sliceCells(k)).FracATFirstLast.sumActiveTime;
            FracATFirstLastOutput{9,k} = EventData(sliceCells(k)).FracATFirstLast.fracActiveTime;
            if isnan(FracATFirstLastOutput{9,k})
                FracATFirstLastOutput{1,k} = NaN;
                FracATFirstLastOutput{3,k} = NaN;
                FracATFirstLastOutput{5,k} = NaN;
                FracATFirstLastOutput{7,k} = NaN;
                FracATFirstLastOutput{9,k} = NaN;
            else
                TotalCells = TotalCells + 1;
            end
        end  
    end

    
    % Calculate means
    FracATFirstLastHeaders{1,end-1} = 'Mean';
    nQuiescentPeriodsMean = nanmean(cell2mat(FracATFirstLastOutput(1,:))); FracATFirstLastOutput{1,end-1} = nQuiescentPeriodsMean;
    sumQuiescentTimeMean  = nanmean(cell2mat(FracATFirstLastOutput(3,:))); FracATFirstLastOutput{3,end-1} = sumQuiescentTimeMean;
    fracQuiescentTimeMean = nanmean(cell2mat(FracATFirstLastOutput(5,:))); FracATFirstLastOutput{5,end-1} = fracQuiescentTimeMean;
    sumActiveTimeMean     = nanmean(cell2mat(FracATFirstLastOutput(7,:))); FracATFirstLastOutput{7,end-1} = sumActiveTimeMean;
    fracActiveTimeMean    = nanmean(cell2mat(FracATFirstLastOutput(9,:))); FracATFirstLastOutput{9,end-1} = fracActiveTimeMean;
    
    % Calculate SE
    FracATFirstLastHeaders{1,end} = 'SE';
    nQuiescentPeriodsSE   = nanstd(cell2mat(FracATFirstLastOutput(1,:)))./sqrt(TotalCells); FracATFirstLastOutput{1,end} = nQuiescentPeriodsSE;
    sumQuiescentTimeSE    = nanstd(cell2mat(FracATFirstLastOutput(3,:)))./sqrt(TotalCells); FracATFirstLastOutput{3,end} = sumQuiescentTimeSE;
    fracQuiescentTimeSE   = nanstd(cell2mat(FracATFirstLastOutput(5,:)))./sqrt(TotalCells); FracATFirstLastOutput{5,end} = fracQuiescentTimeSE;
    sumActiveTimeSE       = nanstd(cell2mat(FracATFirstLastOutput(7,:)))./sqrt(TotalCells); FracATFirstLastOutput{7,end} = sumActiveTimeSE;
    fracActiveTimeSE      = nanstd(cell2mat(FracATFirstLastOutput(9,:)))./sqrt(TotalCells); FracATFirstLastOutput{9,end} = fracActiveTimeSE;
    
    FracATFirstLastOutput = [FracATFirstLastHeaders; FracATFirstLastOutput];
    FracATFirstLastLabels = cell(10,1);
    FracATFirstLastLabels{1,1}  = 'Cell #';
    FracATFirstLastLabels{2,1}  = '# Quiescent Periods';
    FracATFirstLastLabels{4,1}  = 'Sum Quiescent Time';
    FracATFirstLastLabels{6,1}  = 'Fractional Quiescent Time';
    FracATFirstLastLabels{8,1}  = 'Sum Active Time';
    FracATFirstLastLabels{10,1} = 'Fractional Active Time';
    
    FracATFirstLastOutput = [FracATFirstLastLabels, FracATFirstLastOutput];
    
    if ispc
        xlswrite(sliceFile, FracATFirstLastOutput, 'FracActiveFirstLast');
    else
        xlwrite(sliceFile, FracATFirstLastOutput, 'FracActiveFirstLast');
    end
end