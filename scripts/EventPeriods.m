function EventPeriodsOutput = EventPeriods(currCell, threshold)


    if ~isfield(currCell, 'AllPeriods')
        currCell.AllPeriods = AllPeriods(currCell); 
    end

    if threshold == 0
        EventPeriodsOutput = currCell.AllPeriods;
    else
        currCell.AllPeriods = currCell.AllPeriods;
        EventPeriodsOutput = currCell.AllPeriods(currCell.AllPeriods <= threshold);
    end


end