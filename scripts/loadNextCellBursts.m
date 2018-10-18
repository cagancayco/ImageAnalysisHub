function loadNextCellBursts(app)
    load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));
    ax = app.UIAxes;
    cellIdx = app.CellCounter;
    title(ax,app.CellDropDown.Items{cellIdx},'Interpreter','none');
    
    plot(ax,EventData(cellIdx).time, EventData(cellIdx).trace)
    
    if strcmp(EventData(cellIdx).algorithm, 'letsDetect') && ~isempty(EventData(cellIdx).Bursts)
        
        idx = zeros(1, numel(EventData(cellIdx).events));
        for j = 1:length(idx)
            idx(j) = find(EventData(cellIdx).time == EventData(cellIdx).events(j));
        end
        
        [burstID, burstIdx, ~] = unique(EventData(cellIdx).Bursts(:,1));
        
        burstXmin = zeros(numel(burstID),1);
        burstXmax = zeros(numel(burstID),1);
        burstYmin = zeros(numel(burstID),1);
        burstYmax = zeros(numel(burstID),1);
        
        traceRange = range(EventData(cellIdx).trace);
        
        hold(ax, 'on')
        plot(ax, EventData(cellIdx).events, EventData(cellIdx).trace(idx), 'r*')
        
        for i = 1:length(burstID)
            burstXmin(i) = EventData(cellIdx).Bursts(burstIdx(i),2);
            
            if i ~= numel(burstID)
                burstXmax(i) = EventData(cellIdx).Bursts(burstIdx(i+1)-1,2);
            else
                burstXmax(i) = EventData(cellIdx).Bursts(end,2);
            end
            
            amps = EventData(cellIdx).trace(find(EventData(cellIdx).time == burstXmin(i)):find(EventData(cellIdx).time == burstXmax(i)));
            
            burstYmin(i) = min(amps);
            burstYmax(i) = max(amps);
            
            r1 = burstXmin(i)-2;
            r2 = burstYmin(i) - traceRange * 0.05;
            r3 = burstXmax(i) + 2 - r1;
            r4 = burstYmax(i) + traceRange * 0.05 - r2;
            
            rectangle(ax, 'Position', [r1 r2 r3 r4], 'Curvature', 0.2, 'EdgeColor', 'k', 'LineWidth', 2)
            
        end

        hold(ax, 'off')
        
    elseif strcmp(EventData(cellIdx).algorithm, 'minEASE') && ~isempty(EventData(cellIdx).Bursts)
        
        eventsToInc = zeros(length(EventData(cellIdx).events),1);
        for k = 1:length(eventsToInc)
            if ~isempty(find(EventData(cellIdx).classesToInc == EventData(cellIdx).eventClass(k),1))
                eventsToInc(k) = 1;
            end
        end
        eventsToInc = logical(eventsToInc);
        
        [burstID, burstIdx, ~] = unique(EventData(cellIdx).Bursts(:,1));
        
        burstXmin = zeros(numel(burstID),1);
        burstXmax = zeros(numel(burstID),1);
        burstYmin = zeros(numel(burstID),1);
        burstYmax = zeros(numel(burstID),1);
        
        traceRange = range(EventData(cellIdx).trace);
        
        hold(ax, 'on')
        plot(ax, EventData(cellIdx).time(EventData(cellIdx).events(eventsToInc)), EventData(cellIdx).trace(EventData(cellIdx).events(eventsToInc)),'r*')
        for i = 1:length(burstID)
            burstXmin(i) = EventData(cellIdx).Bursts(burstIdx(i),2);
            
            if i ~= numel(burstID)
                burstXmax(i) = EventData(cellIdx).Bursts(burstIdx(i+1)-1,2);
            else
                burstXmax(i) = EventData(cellIdx).Bursts(end,2);
            end
            
            
            amps = EventData(cellIdx).trace(find(EventData(cellIdx).time == burstXmin(i)):find(EventData(cellIdx).time == burstXmax(i)));
            
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
    ax.XLim = [EventData(cellIdx).startTime, EventData(cellIdx).endTime];            
    ax.YLimMode = 'auto';
end