function FracATWindowsOutput = FracATWindows(currCell, threshold)

    
    if ~isfield(currCell, 'AllPeriods')
        currCell.AllPeriods = AllPeriods(currCell); 
    end

    if ~isfield(currCell, 'EventPeriods')
        currCell.EventPeriods = EventPeriods(currCell,threshold);
    end

    if ~isfield(currCell, 'QuiescentPeriods')
        currCell.QuiescentPeriods = QuiescentPeriods(currCell,threshold);
    end

    FracATWindowsOutput = struct;

    FracATWindowsOutput.nQuiescentPeriods = numel(currCell.QuiescentPeriods);
    FracATWindowsOutput.sumQuiescentTime  = sum([currCell.QuiescentPeriods]);
    FracATWindowsOutput.fracQuiescentTime = sum([currCell.QuiescentPeriods])/sum([currCell.AllPeriods]);

    FracATWindowsOutput.sumActiveTime     = sum([currCell.EventPeriods]);
    FracATWindowsOutput.fracActiveTime    = sum([currCell.EventPeriods])/sum([currCell.AllPeriods]);
    
    if FracATWindowsOutput.sumQuiescentTime == 0 && FracATWindowsOutput.sumActiveTime == 0
        FracATWindowsOutput.nQuiescentPeriods = 1;
        FracATWindowsOutput.sumQuiescentTime  = currCell.endTime - currCell.startTime;
        FracATWindowsOutput.fracQuiescentTime = 1;
        FracATWindowsOutput.fracActiveTime    = 0;
    end
        
    if ~isempty(currCell.analysisWindows)
        if ~isfield(currCell, 'AllPeriods_Windows')
            currCell.AllPeriods_Windows = AllPeriods_Windows(currCell);
        end
        
        if ~isfield(currCell, 'EventPeriods_Windows')
            currCell.EventPeriods_Windows = EventPeriods_Windows(currCell,threshold);
        end

        if ~isfield(currCell, 'QuiescentPeriods_Windows')
            currCell.QuiescentPeriods_Windows = QuiescentPeriods_Windows(currCell,threshold);
        end
        
        FracATWindowsOutput.Window1.nQuiescentPeriods = numel(currCell.QuiescentPeriods_Windows.Window1);
        FracATWindowsOutput.Window1.sumQuiescentTime  = sum([currCell.QuiescentPeriods_Windows.Window1]);
        FracATWindowsOutput.Window1.fracQuiescentTime = sum([currCell.QuiescentPeriods_Windows.Window1])/sum([currCell.AllPeriods_Windows.Window1]);
    
        FracATWindowsOutput.Window1.sumActiveTime     = sum([currCell.EventPeriods_Windows.Window1]);
        FracATWindowsOutput.Window1.fracActiveTime    = sum([currCell.EventPeriods_Windows.Window1])/sum([currCell.AllPeriods_Windows.Window1]);    
        
        
        FracATWindowsOutput.Window2.nQuiescentPeriods = numel(currCell.QuiescentPeriods_Windows.Window2);
        FracATWindowsOutput.Window2.sumQuiescentTime  = sum([currCell.QuiescentPeriods_Windows.Window2]);
        FracATWindowsOutput.Window2.fracQuiescentTime = sum([currCell.QuiescentPeriods_Windows.Window2])/sum([currCell.AllPeriods_Windows.Window2]);
    
        FracATWindowsOutput.Window2.sumActiveTime     = sum([currCell.EventPeriods_Windows.Window2]);
        FracATWindowsOutput.Window2.fracActiveTime    = sum([currCell.EventPeriods_Windows.Window2])/sum([currCell.AllPeriods_Windows.Window2]);  
    end


end