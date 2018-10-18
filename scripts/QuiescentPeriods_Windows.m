function QuiescentPeriods_WindowsOutput = QuiescentPeriods_Windows(currCell, threshold)

    if ~isfield(currCell, 'AllPeriods_Windows')
        currCell.AllPeriods_Windows = AllPeriods_Windows(currCell);
    end
    
    QuiescentPeriods_WindowsOutput = struct;
    
    if threshold == 0
        QuiescentPeriods_WindowsOutput.Window1 = [];
        QuiescentPeriods_WindowsOutput.Window2 = [];
    else
        QuiescentPeriods_WindowsOutput.Window1 = currCell.AllPeriods_Windows.Window1(currCell.AllPeriods_Windows.Window1 > threshold);
        QuiescentPeriods_WindowsOutput.Window2 = currCell.AllPeriods_Windows.Window2(currCell.AllPeriods_Windows.Window2 > threshold);
    end


end