function AllPeriods_WindowsOutput = AllPeriods_Windows(currCell)

    if strcmp(currCell.algorithm,'minEASE')
        events = currCell.events;
        events = currCell.time(events)';
        
        classesToInc = currCell.classesToInc;
        eventClass   = currCell.eventClass;
        eventsToInc = zeros(length(events),1);
        for i = 1:length(eventsToInc)
            if ~isempty(find(classesToInc == eventClass(i),1))
                eventsToInc(i) = 1;
            end
        end
        eventsToInc = logical(eventsToInc);
        
        events = events(eventsToInc);
    else
        events = currCell.events;
    end

    win1_MinTime = currCell.analysisWindows.Var2(1);
    win1_MaxTime = currCell.analysisWindows.Var3(1);


    win1_Events  = events;
    win1_Events(win1_Events < win1_MinTime) = [];
    win1_Events(win1_Events > win1_MaxTime) = [];

    win1_AllPeriods  = diff(win1_Events);
    if ~isempty(win1_AllPeriods)
        win1_firstPeriod = win1_Events(1) - win1_MinTime;
        win1_lastPeriod  = win1_MaxTime - win1_Events(end);
        
        if win1_firstPeriod ~= 0
            win1_AllPeriods = [win1_firstPeriod; win1_AllPeriods];
        end
        
        if win1_lastPeriods ~= 0
            win1_AllPeriods = [win1_AllPeriods; win1_lastPeriod];
        end
    end

    win2_MinTime = currCell.analysisWindows.Var2(2);
    win2_MaxTime = currCell.analysisWindows.Var3(2);


    win2_Events  = events;
    win2_Events(win2_Events < win2_MinTime) = [];
    win2_Events(win2_Events > win2_MaxTime) = [];

    win2_AllPeriods  = diff(win2_Events);
    if ~isempty(win2_AllPeriods)
        win2_firstPeriod = win2_Events(1) - win2_MinTime;
        win2_lastPeriod  = win2_MaxTime - win2_Events(end);
        
        if win2_firstPeriod ~= 0
            win2_AllPeriods = [win2_firstPeriod; win2_AllPeriods];
        end
        
        if win2_lastPeriods ~= 0
            win2_AllPeriods = [win2_AllPeriods; win2_lastPeriod];
        end
    end

    AllPeriods_WindowsOutput = struct;
    AllPeriods_WindowsOutput.Window1 = win1_AllPeriods;
    AllPeriods_WindowsOutput.Window2 = win2_AllPeriods;
end