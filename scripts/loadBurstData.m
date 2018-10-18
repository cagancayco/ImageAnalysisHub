function loadBurstData(app)
    load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value))
    
    
    for i = 1:length(EventData)
        app.CellDropDown.Items{i} = EventData(i).cellName;
        %EventData(i).Bursts = FindBursts(EventData(i));
    end
    
    ax = app.UIAxes;
    title(ax,app.CellDropDown.Value,'Interpreter','none');
    plot(ax, EventData(1).time, EventData(1).trace)
    
    if strcmp(EventData(1).algorithm, 'letsDetect')
        
        idx = zeros(1, numel(EventData(1).events));
        for j = 1:length(idx)
            idx(j) = find(EventData(1).time == EventData(1).events(j));
        end
        
        [burstID, burstIdx, ~] = unique(EventData(1).Bursts(:,1));
        
        burstXmin = zeros(numel(burstID),1);
        burstXmax = zeros(numel(burstID),1);
        burstYmin = zeros(numel(burstID),1);
        burstYmax = zeros(numel(burstID),1);
        
        traceRange = range(EventData(1).trace);
        
        hold(ax, 'on')
        plot(ax, EventData(1).events, EventData(1).trace(idx), 'r*')
        
        for i = 1:length(burstID)
            burstXmin(i) = EventData(1).Bursts(burstIdx(i),2);
            
            if i ~= numel(burstID)
                burstXmax(i) = EventData(1).Bursts(burstIdx(i+1)-1,2);
            else
                burstXmax(i) = EventData(1).Bursts(end,2);
            end
            
            amps = EventData(1).trace(find(EventData(1).time == burstXmin(i)):find(EventData(1).time == burstXmax(i)));
            
            burstYmin(i) = min(amps);
            burstYmax(i) = max(amps);
            
            r1 = burstXmin(i)-2;
            r2 = burstYmin(i) - traceRange * 0.05;
            r3 = burstXmax(i) + 2 - r1;
            r4 = burstYmax(i) + traceRange * 0.05 - r2;
            
            rectangle(ax, 'Position', [r1 r2 r3 r4], 'Curvature', 0.2, 'EdgeColor', 'k', 'LineWidth', 2)
            
        end

        hold(ax, 'off')
        
    else
        
        eventsToInc = zeros(length(EventData(1).events),1);
        for k = 1:length(eventsToInc)
            if ~isempty(find(EventData(1).classesToInc == EventData(1).eventClass(k),1))
                eventsToInc(k) = 1;
            end
        end
        eventsToInc = logical(eventsToInc);
        
        [burstID, burstIdx, ~] = unique(EventData(1).Bursts(:,1));
        
        burstXmin = zeros(numel(burstID),1);
        burstXmax = zeros(numel(burstID),1);
        burstYmin = zeros(numel(burstID),1);
        burstYmax = zeros(numel(burstID),1);
        
        traceRange = range(EventData(1).trace);
        
        hold(ax, 'on')
        plot(ax, EventData(1).time(EventData(1).events(eventsToInc)), EventData(1).trace(EventData(1).events(eventsToInc)),'r*')
        for i = 1:length(burstID)
            burstXmin(i) = EventData(1).Bursts(burstIdx(i),2);
            
            if i ~= numel(burstID)
                burstXmax(i) = EventData(1).Bursts(burstIdx(i+1)-1,2);
            else
                burstXmax(i) = EventData(1).Bursts(end,2);
            end
            
            
            amps = EventData(1).trace(find(EventData(1).time == burstXmin(i)):find(EventData(1).time == burstXmax(i)));
            
            burstYmin(i) = min(amps);
            burstYmax(i) = max(amps);
            
            r1 = burstXmin(i)-1;
            r2 = burstYmin(i) - traceRange * 0.05;
            r3 = burstXmax(i) + 1 - r1;
            r4 = burstYmax(i) + traceRange * 0.05 - r2;
            
            rectangle(ax, 'Position', [r1 r2 r3 r4], 'Curvature', 0.2, 'EdgeColor', 'k', 'LineWidth', 2)
            
        end
        hold(ax, 'off')
        
    end
    ax.XLim = [EventData(1).startTime, EventData(1).endTime];            
    ax.YLimMode = 'auto';
end