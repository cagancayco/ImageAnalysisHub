function loadNextCell(app)

    load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));
    ax = app.UIAxes;
    cellIdx = app.CellCounter;
    title(ax,app.CellDropDown.Items{cellIdx},'Interpreter','none');
    
    plot(ax,EventData(cellIdx).time, EventData(cellIdx).trace)
    if strcmp(EventData(cellIdx).algorithm, 'letsDetect')
        

        idx = zeros(1,numel(EventData(cellIdx).events));
        for j = 1:length(idx)
            idx(j) = find(EventData(cellIdx).time == EventData(cellIdx).events(j));
        end
        hold(ax, 'on')
        plot(ax,EventData(cellIdx).events, EventData(cellIdx).trace(idx), 'r*')
        hold(ax, 'off')
        
    else
        
        eventsToInc = zeros(length(EventData(cellIdx).events),1);
        for k = 1:length(eventsToInc)
            if ~isempty(find(EventData(cellIdx).classesToInc == EventData(cellIdx).eventClass(k),1))
                eventsToInc(k) = 1;
            end
        end
        eventsToInc = logical(eventsToInc);
        
        hold(ax, 'on')
        plot(ax, EventData(cellIdx).time(EventData(cellIdx).events(eventsToInc)), EventData(cellIdx).trace(EventData(cellIdx).events(eventsToInc)),'r*')
        hold(ax, 'off')
        
    end


    ax.XLim = [EventData(cellIdx).startTime, EventData(cellIdx).endTime];            
    ax.YLimMode = 'auto';





end