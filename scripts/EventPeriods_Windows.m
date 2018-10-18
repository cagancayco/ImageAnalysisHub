function EventPeriods_WindowsOutput = EventPeriods_Windows(currCell, threshold)

    if ~isfield(currCell, 'AllPeriods_Windows')
        currCell.AllPeriods_Windows = AllPeriods_Windows(currCell);
    end
    
    EventPeriods_WindowsOutput = struct;
    
    if threshold == 0
        EventPeriods_WindowsOutput.Window1 = currCell.AllPeriods_Windows.Window1;
        EventPeriods_WindowsOutput.Window2 = currCell.AllPeriods_Windows.Window2;
    else
        EventPeriods_WindowsOutput.Window1 = currCell.AllPeriods_Windows.Window1(currCell.AllPeriods_Windows.Window1 <= threshold);
        EventPeriods_WindowsOutput.Window2 = currCell.AllPeriods_Windows.Window2(currCell.AllPeriods_Windows.Window2 <= threshold);
    end


end