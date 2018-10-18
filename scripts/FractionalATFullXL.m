function FractionalATFullXL(EventData, currSlice, sliceCells, sliceFile)

    FracATWindowsFullHeaders = cell(1,numel(sliceCells)+2);
    FracATWindowsFullOutput  = cell(9,numel(sliceCells)+2);
    TotalCells = 0;

    for k = 1:numel(sliceCells)
        FracATWindowsFullHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).FracATWindows)
            FracATWindowsFullOutput{1,k} = EventData(sliceCells(k)).FracATWindows.nQuiescentPeriods;
            FracATWindowsFullOutput{3,k} = EventData(sliceCells(k)).FracATWindows.sumQuiescentTime;
            FracATWindowsFullOutput{5,k} = EventData(sliceCells(k)).FracATWindows.fracQuiescentTime;
            FracATWindowsFullOutput{7,k} = EventData(sliceCells(k)).FracATWindows.sumActiveTime;
            FracATWindowsFullOutput{9,k} = EventData(sliceCells(k)).FracATWindows.fracActiveTime;
            TotalCells = TotalCells + 1;
        end  
    end

    
    % Calculate means
    FracATWindowsFullHeaders{1,end-1} = 'Mean';
    nQuiescentPeriodsMean = nanmean(cell2mat(FracATWindowsFullOutput(1,:))); FracATWindowsFullOutput{1,end-1} = nQuiescentPeriodsMean;
    sumQuiescentTimeMean  = nanmean(cell2mat(FracATWindowsFullOutput(3,:))); FracATWindowsFullOutput{3,end-1} = sumQuiescentTimeMean;
    fracQuiescentTimeMean = nanmean(cell2mat(FracATWindowsFullOutput(5,:))); FracATWindowsFullOutput{5,end-1} = fracQuiescentTimeMean;
    sumActiveTimeMean     = nanmean(cell2mat(FracATWindowsFullOutput(7,:))); FracATWindowsFullOutput{7,end-1} = sumActiveTimeMean;
    fracActiveTimeMean    = nanmean(cell2mat(FracATWindowsFullOutput(9,:))); FracATWindowsFullOutput{9,end-1} = fracActiveTimeMean;
    
    % Calculate SE
    FracATWindowsFullHeaders{1,end} = 'SE';
    nQuiescentPeriodsSE   = nanstd(cell2mat(FracATWindowsFullOutput(1,:)))./sqrt(TotalCells); FracATWindowsFullOutput{1,end} = nQuiescentPeriodsSE;
    sumQuiescentTimeSE    = nanstd(cell2mat(FracATWindowsFullOutput(3,:)))./sqrt(TotalCells); FracATWindowsFullOutput{3,end} = sumQuiescentTimeSE;
    fracQuiescentTimeSE   = nanstd(cell2mat(FracATWindowsFullOutput(5,:)))./sqrt(TotalCells); FracATWindowsFullOutput{5,end} = fracQuiescentTimeSE;
    sumActiveTimeSE       = nanstd(cell2mat(FracATWindowsFullOutput(7,:)))./sqrt(TotalCells); FracATWindowsFullOutput{7,end} = sumActiveTimeSE;
    fracActiveTimeSE      = nanstd(cell2mat(FracATWindowsFullOutput(9,:)))./sqrt(TotalCells); FracATWindowsFullOutput{9,end} = fracActiveTimeSE;
    
    FracATWindowsFullOutput = [FracATWindowsFullHeaders; FracATWindowsFullOutput];
    FracATWindowsFullLabels = cell(10,1);
    FracATWindowsFullLabels{1,1}  = 'Cell #';
    FracATWindowsFullLabels{2,1}  = '# Quiescent Periods';
    FracATWindowsFullLabels{4,1}  = 'Sum Quiescent Time';
    FracATWindowsFullLabels{6,1}  = 'Fractional Quiescent Time';
    FracATWindowsFullLabels{8,1}  = 'Sum Active Time';
    FracATWindowsFullLabels{10,1} = 'Fractional Active Time';
    
    FracATWindowsFullOutput = [FracATWindowsFullLabels, FracATWindowsFullOutput];
    
    if ispc
        xlswrite(sliceFile, FracATWindowsFullOutput, 'FracActiveFullTime');
    else
        xlwrite(sliceFile, FracATWindowsFullOutput, 'FracActiveFullTime');
    end






end