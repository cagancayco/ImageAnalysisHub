function EventPeriodsXL(EventData, currSlice, sliceCells, sliceFile)


    EventPeriodsHeaders = cell(1,numel(sliceCells));
    EventPeriodsResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    EventPeriodsStats = cell(6, numel(sliceCells) + 2);
    EventPeriodsStats{1,1} = 'mean';
    EventPeriodsStats{2,1} = 'SE';
    EventPeriodsStats{3,1} = 'median';
    EventPeriodsStats{4,1} = 'max';
    EventPeriodsStats{5,1} = 'min';
    EventPeriodsStats{6,1} = 'N';
            
    for k = 1:numel(sliceCells)
        EventPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).EventPeriods)
            EventPeriodsResults{k} = [EventData(sliceCells(k)).EventPeriods];
            rowCounts(k)  = numel([EventData(sliceCells(k)).EventPeriods]);
            
            % Mean
            EventPeriodsStats{1,k+1} = mean(EventPeriodsResults{k});
            % SE
            EventPeriodsStats{2,k+1} = std(EventPeriodsResults{k})./sqrt(rowCounts(k));
            % Median
            EventPeriodsStats{3,k+1} = median(EventPeriodsResults{k});
            % Max
            EventPeriodsStats{4,k+1} = max(EventPeriodsResults{k});
            % Min
            EventPeriodsStats{5,k+1} = min(EventPeriodsResults{k});
            % N
            EventPeriodsStats{6,k+1} = rowCounts(k);
        else
            rowCounts(k) = 0;
        end
    end
            
    EventPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        EventPeriodsOutput(1:rowCounts(k),k) = num2cell(EventPeriodsResults{k});
    end
    TotalPeriods = reshape(EventPeriodsOutput, [numel(EventPeriodsOutput) ,1]);
    TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
    EventPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));
    
    % Mean
    EventPeriodsStats{1,end} = mean(TotalPeriods);
    % SE
    EventPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
    % Median
    EventPeriodsStats{3,end} = median(TotalPeriods);
    % Max
    EventPeriodsStats{4,end} = max(TotalPeriods);
    % Min
    EventPeriodsStats{5,end} = min(TotalPeriods);
    % N
    EventPeriodsStats{6,end} = numel(TotalPeriods);
            
    EventPeriodsOutput = [[{'Descriptive Stats'}, EventPeriodsHeaders, {'Event Periods'}]; EventPeriodsStats; [{''}, EventPeriodsHeaders, {'Event Periods'}];[cell(numel(TotalPeriods),1), EventPeriodsOutput]];
    
    if ispc
        xlswrite(sliceFile, EventPeriodsOutput, 'Event Periods');
    else
        xlwrite(sliceFile, EventPeriodsOutput, 'Event Periods');
    end        
    
    %% Window Analysis
    if ~isempty(EventData(1).analysisWindows)
        % Window 1
        EventPeriodsHeaders = cell(1,numel(sliceCells));
        EventPeriodsResults = cell(1,numel(sliceCells));
        rowCounts = zeros(1,numel(sliceCells));
    
        % Calculate stats
        EventPeriodsStats = cell(6, numel(sliceCells) + 2);
        EventPeriodsStats{1,1} = 'mean';
        EventPeriodsStats{2,1} = 'SE';
        EventPeriodsStats{3,1} = 'median';
        EventPeriodsStats{4,1} = 'max';
        EventPeriodsStats{5,1} = 'min';
        EventPeriodsStats{6,1} = 'N';
        

        for k = 1:numel(sliceCells)
            EventPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            if ~isempty(EventData(sliceCells(k)).EventPeriods_Windows)
                EventPeriodsResults{k} = [EventData(sliceCells(k)).EventPeriods_Windows.Window1];
                rowCounts(k)  = numel([EventData(sliceCells(k)).EventPeriods_Windows.Window1]);

                % Mean
                EventPeriodsStats{1,k+1} = mean(EventPeriodsResults{k});
                % SE
                EventPeriodsStats{2,k+1} = std(EventPeriodsResults{k})./sqrt(rowCounts(k));
                % Median
                EventPeriodsStats{3,k+1} = median(EventPeriodsResults{k});
                % Max
                EventPeriodsStats{4,k+1} = max(EventPeriodsResults{k});
                % Min
                EventPeriodsStats{5,k+1} = min(EventPeriodsResults{k});
                % N
                EventPeriodsStats{6,k+1} = rowCounts(k);
            else
                rowCounts(k) = 0;
            end
        end

        EventPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
        for k = 1:numel(sliceCells)
            EventPeriodsOutput(1:rowCounts(k),k) = num2cell(EventPeriodsResults{k});
        end
        TotalPeriods = reshape(EventPeriodsOutput, [numel(EventPeriodsOutput) ,1]);
        TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
        EventPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));

        % Mean
        EventPeriodsStats{1,end} = mean(TotalPeriods);
        % SE
        EventPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
        % Median
        EventPeriodsStats{3,end} = median(TotalPeriods);
        % Max
        EventPeriodsStats{4,end} = max(TotalPeriods);
        % Min
        EventPeriodsStats{5,end} = min(TotalPeriods);
        % N
        EventPeriodsStats{6,end} = numel(TotalPeriods);

        EventPeriodsOutput = [[{'Descriptive Stats'}, EventPeriodsHeaders, {'Event Periods Window1'}]; EventPeriodsStats; [{''}, EventPeriodsHeaders, {'Event Periods Window1'}];[cell(numel(TotalPeriods),1), EventPeriodsOutput]];
        
        if ispc
            xlswrite(sliceFile, EventPeriodsOutput, 'Event Periods Window1');
        else
            xlwrite(sliceFile, EventPeriodsOutput, 'Event Periods Window1');
        end
        
        % Window 2
        EventPeriodsHeaders = cell(1,numel(sliceCells));
        EventPeriodsResults = cell(1,numel(sliceCells));
        rowCounts = zeros(1,numel(sliceCells));
    
        % Calculate stats
        EventPeriodsStats = cell(6, numel(sliceCells) + 2);
        EventPeriodsStats{1,1} = 'mean';
        EventPeriodsStats{2,1} = 'SE';
        EventPeriodsStats{3,1} = 'median';
        EventPeriodsStats{4,1} = 'max';
        EventPeriodsStats{5,1} = 'min';
        EventPeriodsStats{6,1} = 'N';

        for k = 1:numel(sliceCells)
            EventPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            if ~isempty(EventData(sliceCells(k)).EventPeriods_Windows)
                EventPeriodsResults{k} = [EventData(sliceCells(k)).EventPeriods_Windows.Window2];
                rowCounts(k)  = numel([EventData(sliceCells(k)).EventPeriods_Windows.Window2]);

                % Mean
                EventPeriodsStats{1,k+1} = mean(EventPeriodsResults{k});
                % SE
                EventPeriodsStats{2,k+1} = std(EventPeriodsResults{k})./sqrt(rowCounts(k));
                % Median
                EventPeriodsStats{3,k+1} = median(EventPeriodsResults{k});
                % Max
                EventPeriodsStats{4,k+1} = max(EventPeriodsResults{k});
                % Min
                EventPeriodsStats{5,k+1} = min(EventPeriodsResults{k});
                % N
                EventPeriodsStats{6,k+1} = rowCounts(k);
            else
                rowCounts(k) = 0;
            end
        end

        EventPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
        for k = 1:numel(sliceCells)
            EventPeriodsOutput(1:rowCounts(k),k) = num2cell(EventPeriodsResults{k});
        end
        TotalPeriods = reshape(EventPeriodsOutput, [numel(EventPeriodsOutput) ,1]);
        TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
        EventPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));

        % Mean
        EventPeriodsStats{1,end} = mean(TotalPeriods);
        % SE
        EventPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
        % Median
        EventPeriodsStats{3,end} = median(TotalPeriods);
        % Max
        EventPeriodsStats{4,end} = max(TotalPeriods);
        % Min
        EventPeriodsStats{5,end} = min(TotalPeriods);
        % N
        EventPeriodsStats{6,end} = numel(TotalPeriods);

        EventPeriodsOutput = [[{'Descriptive Stats'}, EventPeriodsHeaders, {'Event Periods Window2'}]; EventPeriodsStats; [{''}, EventPeriodsHeaders, {'Event Periods Window2'}];[cell(numel(TotalPeriods),1), EventPeriodsOutput]];
        
        if ispc
            xlswrite(sliceFile, EventPeriodsOutput, 'Event Periods Window2');
        else
            xlwrite(sliceFile, EventPeriodsOutput, 'Event Periods Window2');
        end
    end










end