function loadTraceData(app)
    load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));       
    
    
    for i = 1:length(EventData)
        app.CellDropDown.Items{i} = EventData(i).cellName;
    end

    
    
    ax = app.UIAxes;
    title(ax,app.CellDropDown.Value,'Interpreter','none');
    plot(ax, EventData(1).time, EventData(1).trace)
    
    if strcmp(EventData(1).algorithm, 'letsDetect')
        
                    
        idx = zeros(1,numel(EventData(1).events));
        for j = 1:length(idx)
            idx(j) = find(EventData(1).time == EventData(1).events(j));
        end
        hold(ax, 'on')
        plot(ax,EventData(1).events, EventData(1).trace(idx), 'r*')
        hold(ax, 'off')
        
    else
        
        eventsToInc = zeros(length(EventData(1).events),1);
        for k = 1:length(eventsToInc)
            if ~isempty(find(EventData(1).classesToInc == EventData(1).eventClass(k),1))
                eventsToInc(k) = 1;
            end
        end
        eventsToInc = logical(eventsToInc);
        
        hold(ax, 'on')
        plot(ax, EventData(1).time(EventData(1).events(eventsToInc)), EventData(1).trace(EventData(1).events(eventsToInc)),'r*')
        hold(ax, 'off')
                   
    end


    ax.XLim = [EventData(1).startTime, EventData(1).endTime];            
    ax.YLimMode = 'auto'; 


end