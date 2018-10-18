function FracATFirstLastOutput = FracATFirstLast(currCell, threshold)

    if ~isfield(currCell, 'AllPeriods')
        currCell.AllPeriods = AllPeriods(currCell); 
    end
    
    currCell.AllPeriods = currCell.AllPeriods(2:end-1);
    
    if threshold == 0
        currCell.EventPeriods = currCell.AllPeriods;
    else
        currCell.EventPeriods = currCell.AllPeriods(currCell.AllPeriods <= threshold);
    end

    if threshold == 0
        currCell.QuiescentPeriods = [];
    else
        currCell.QuiescentPeriods = currCell.AllPeriods(currCell.AllPeriods > threshold);
    end
    
    FracATFirstLastOutput = struct;
    
    FracATFirstLastOutput.nQuiescentPeriods = numel(currCell.QuiescentPeriods);
    FracATFirstLastOutput.sumQuiescentTime  = sum([currCell.QuiescentPeriods]);
    FracATFirstLastOutput.fracQuiescentTime = sum([currCell.QuiescentPeriods])/sum([currCell.AllPeriods]);
    
    FracATFirstLastOutput.sumActiveTime     = sum([currCell.EventPeriods]);
    FracATFirstLastOutput.fracActiveTime    = sum([currCell.EventPeriods])/sum([currCell.AllPeriods]);


end