function FractionalATWindowsXL(EventData, currSlice, sliceCells, sliceFile)
    
    FracATWindowsHeaders  = cell(1, numel(sliceCells)+2);
    FracATWindowsOutput = cell(19, numel(sliceCells)+2);
    TotalCells = 0;
    
    for k = 1:numel(sliceCells)
        FracATWindowsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).FracATWindows)
            FracATWindowsOutput{2,k}  = EventData(sliceCells(k)).FracATWindows.Window1.nQuiescentPeriods;
            FracATWindowsOutput{3,k}  = EventData(sliceCells(k)).FracATWindows.Window2.nQuiescentPeriods;
            FracATWindowsOutput{6,k}  = EventData(sliceCells(k)).FracATWindows.Window1.sumQuiescentTime;
            FracATWindowsOutput{7,k}  = EventData(sliceCells(k)).FracATWindows.Window2.sumQuiescentTime;
            FracATWindowsOutput{10,k} = EventData(sliceCells(k)).FracATWindows.Window1.fracQuiescentTime;
            FracATWindowsOutput{11,k} = EventData(sliceCells(k)).FracATWindows.Window2.fracQuiescentTime;
            FracATWindowsOutput{14,k} = EventData(sliceCells(k)).FracATWindows.Window1.sumActiveTime;
            FracATWindowsOutput{15,k} = EventData(sliceCells(k)).FracATWindows.Window2.sumActiveTime;
            FracATWindowsOutput{18,k} = EventData(sliceCells(k)).FracATWindows.Window1.fracActiveTime;
            FracATWindowsOutput{19,k} = EventData(sliceCells(k)).FracATWindows.Window2.fracActiveTime;
            TotalCells = TotalCells + 1;
        end   
    end

    % Calculate means
    FracATWindowsHeaders{1,end-1}  = 'Mean';
    FracATWindowsOutput{2,end-1}   = nanmean(cell2mat(FracATWindowsOutput(2,:)));
    FracATWindowsOutput{3,end-1}   = nanmean(cell2mat(FracATWindowsOutput(3,:)));
    FracATWindowsOutput{6,end-1}   = nanmean(cell2mat(FracATWindowsOutput(6,:)));
    FracATWindowsOutput{7,end-1}   = nanmean(cell2mat(FracATWindowsOutput(7,:)));
    FracATWindowsOutput{10,end-1}  = nanmean(cell2mat(FracATWindowsOutput(10,:)));
    FracATWindowsOutput{11,end-1}  = nanmean(cell2mat(FracATWindowsOutput(11,:)));
    FracATWindowsOutput{14,end-1}  = nanmean(cell2mat(FracATWindowsOutput(14,:)));
    FracATWindowsOutput{15,end-1}  = nanmean(cell2mat(FracATWindowsOutput(15,:)));
    FracATWindowsOutput{18,end-1}  = nanmean(cell2mat(FracATWindowsOutput(18,:)));
    FracATWindowsOutput{19,end-1}  = nanmean(cell2mat(FracATWindowsOutput(19,:)));
    
    % Calculate SE
    FracATWindowsHeaders{1,end}  = 'SE';
    FracATWindowsOutput{2,end}   = nanstd(cell2mat(FracATWindowsOutput(2,:)))./sqrt(TotalCells);
    FracATWindowsOutput{3,end}   = nanstd(cell2mat(FracATWindowsOutput(3,:)))./sqrt(TotalCells);
    FracATWindowsOutput{6,end}   = nanstd(cell2mat(FracATWindowsOutput(6,:)))./sqrt(TotalCells);
    FracATWindowsOutput{7,end}   = nanstd(cell2mat(FracATWindowsOutput(7,:)))./sqrt(TotalCells);
    FracATWindowsOutput{10,end}  = nanstd(cell2mat(FracATWindowsOutput(10,:)))./sqrt(TotalCells);
    FracATWindowsOutput{11,end}  = nanstd(cell2mat(FracATWindowsOutput(11,:)))./sqrt(TotalCells);
    FracATWindowsOutput{14,end}  = nanstd(cell2mat(FracATWindowsOutput(14,:)))./sqrt(TotalCells);
    FracATWindowsOutput{15,end}  = nanstd(cell2mat(FracATWindowsOutput(15,:)))./sqrt(TotalCells);
    FracATWindowsOutput{18,end}  = nanstd(cell2mat(FracATWindowsOutput(18,:)))./sqrt(TotalCells);
    FracATWindowsOutput{19,end}  = nanstd(cell2mat(FracATWindowsOutput(19,:)))./sqrt(TotalCells);
    
    FracATWindowsOutput = [FracATWindowsHeaders; FracATWindowsOutput];
    FracATWindowsLabels = cell(17,1);
    FracATWindowsLabels{1,1} = 'Cell #';
    FracATWindowsLabels{2,1} = '# Quiescent Periods';
    FracATWindowsLabels{3,1} = [num2str(EventData(1).analysisWindows.Var2(1)),'-',num2str(EventData(1).analysisWindows.Var3(1))];
    FracATWindowsLabels{4,1} = [num2str(EventData(1).analysisWindows.Var2(2)),'-',num2str(EventData(1).analysisWindows.Var3(2))];
    FracATWindowsLabels{6,1} = 'Sum Quiescent Time';
    FracATWindowsLabels{7,1} = [num2str(EventData(1).analysisWindows.Var2(1)),'-',num2str(EventData(1).analysisWindows.Var3(1))];
    FracATWindowsLabels{8,1} = [num2str(EventData(1).analysisWindows.Var2(2)),'-',num2str(EventData(1).analysisWindows.Var3(2))];
    FracATWindowsLabels{10,1} = 'Fractional Quiescent Time';
    FracATWindowsLabels{11,1} = [num2str(EventData(1).analysisWindows.Var2(1)),'-',num2str(EventData(1).analysisWindows.Var3(1))];
    FracATWindowsLabels{12,1} = [num2str(EventData(1).analysisWindows.Var2(2)),'-',num2str(EventData(1).analysisWindows.Var3(2))];
    FracATWindowsLabels{14,1} = 'Sum Active Time';
    FracATWindowsLabels{15,1} = [num2str(EventData(1).analysisWindows.Var2(1)),'-',num2str(EventData(1).analysisWindows.Var3(1))];
    FracATWindowsLabels{16,1} = [num2str(EventData(1).analysisWindows.Var2(2)),'-',num2str(EventData(1).analysisWindows.Var3(2))];
    FracATWindowsLabels{18,1} = 'Fractional Active Time';
    FracATWindowsLabels{19,1} = [num2str(EventData(1).analysisWindows.Var2(1)),'-',num2str(EventData(1).analysisWindows.Var3(1))];
    FracATWindowsLabels{20,1} = [num2str(EventData(1).analysisWindows.Var2(2)),'-',num2str(EventData(1).analysisWindows.Var3(2))];
    FracATWindowsOutput = [FracATWindowsLabels, FracATWindowsOutput];
    
    if ispc
        xlswrite(sliceFile, FracATWindowsOutput, 'FracActiveWindows');
    else
        xlwrite(sliceFile, FracATWindowsOutput, 'FracActiveWindows');
    end
end