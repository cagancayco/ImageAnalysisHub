function QuiescentPeriodsOutput = QuiescentPeriods(currCell, threshold)


    if ~isfield(currCell, 'AllPeriods')
        currCell.AllPeriods = AllPeriods(currCell); 
    end

    if threshold == 0
        QuiescentPeriodsOutput = [];
    else
        QuiescentPeriodsOutput = currCell.AllPeriods(currCell.AllPeriods > threshold);
    end


end