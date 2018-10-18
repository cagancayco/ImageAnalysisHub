function QuiescentPeriodsXL(EventData, currSlice, sliceCells, sliceFile)


    QuiescentPeriodsHeaders = cell(1,numel(sliceCells));
    QuiescentPeriodsResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    QuiescentPeriodsStats = cell(6, numel(sliceCells) + 2);
    QuiescentPeriodsStats{1,1} = 'mean';
    QuiescentPeriodsStats{2,1} = 'SE';
    QuiescentPeriodsStats{3,1} = 'median';
    QuiescentPeriodsStats{4,1} = 'max';
    QuiescentPeriodsStats{5,1} = 'min';
    QuiescentPeriodsStats{6,1} = 'N';
            
    for k = 1:numel(sliceCells)
        QuiescentPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).QuiescentPeriods)
            QuiescentPeriodsResults{k} = [EventData(sliceCells(k)).QuiescentPeriods];
            rowCounts(k)  = numel([EventData(sliceCells(k)).QuiescentPeriods]);
            
            % Mean
            QuiescentPeriodsStats{1,k+1} = mean(QuiescentPeriodsResults{k});
            % SE
            QuiescentPeriodsStats{2,k+1} = std(QuiescentPeriodsResults{k})./sqrt(rowCounts(k));
            % Median
            QuiescentPeriodsStats{3,k+1} = median(QuiescentPeriodsResults{k});
            % Max
            QuiescentPeriodsStats{4,k+1} = max(QuiescentPeriodsResults{k});
            % Min
            QuiescentPeriodsStats{5,k+1} = min(QuiescentPeriodsResults{k});
            % N
            QuiescentPeriodsStats{6,k+1} = rowCounts(k);
        else
            rowCounts(k) = 0;
        end
    end
            
    QuiescentPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        QuiescentPeriodsOutput(1:rowCounts(k),k) = num2cell(QuiescentPeriodsResults{k});
    end
    TotalPeriods = reshape(QuiescentPeriodsOutput, [numel(QuiescentPeriodsOutput) ,1]);
    TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
    QuiescentPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));
    
    % Mean
    QuiescentPeriodsStats{1,end} = mean(TotalPeriods);
    % SE
    QuiescentPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
    % Median
    QuiescentPeriodsStats{3,end} = median(TotalPeriods);
    % Max
    QuiescentPeriodsStats{4,end} = max(TotalPeriods);
    % Min
    QuiescentPeriodsStats{5,end} = min(TotalPeriods);
    % N
    QuiescentPeriodsStats{6,end} = numel(TotalPeriods);
            
    QuiescentPeriodsOutput = [[{'Descriptive Stats'}, QuiescentPeriodsHeaders, {'Quiescent Periods'}]; QuiescentPeriodsStats; [{''}, QuiescentPeriodsHeaders, {'Quiescent Periods'}];[cell(numel(TotalPeriods),1), QuiescentPeriodsOutput]];
    if ispc
        xlswrite(sliceFile, QuiescentPeriodsOutput, 'Quiescent Periods');
    else
        xlwrite(sliceFile, QuiescentPeriodsOutput, 'Quiescent Periods');
    end    
            
    if ~isempty(EventData(1).analysisWindows)
        % Window 1
        QuiescentPeriodsHeaders = cell(1,numel(sliceCells));
        QuiescentPeriodsResults = cell(1,numel(sliceCells));
        rowCounts = zeros(1,numel(sliceCells));
    
        % Calculate stats
        QuiescentPeriodsStats = cell(6, numel(sliceCells) + 2);
        QuiescentPeriodsStats{1,1} = 'mean';
        QuiescentPeriodsStats{2,1} = 'SE';
        QuiescentPeriodsStats{3,1} = 'median';
        QuiescentPeriodsStats{4,1} = 'max';
        QuiescentPeriodsStats{5,1} = 'min';
        QuiescentPeriodsStats{6,1} = 'N';
        

        for k = 1:numel(sliceCells)
            QuiescentPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            if ~isempty(EventData(sliceCells(k)).QuiescentPeriods_Windows)
                QuiescentPeriodsResults{k} = [EventData(sliceCells(k)).QuiescentPeriods_Windows.Window1];
                rowCounts(k)  = numel([EventData(sliceCells(k)).QuiescentPeriods_Windows.Window1]);

                % Mean
                QuiescentPeriodsStats{1,k+1} = mean(QuiescentPeriodsResults{k});
                % SE
                QuiescentPeriodsStats{2,k+1} = std(QuiescentPeriodsResults{k})./sqrt(rowCounts(k));
                % Median
                QuiescentPeriodsStats{3,k+1} = median(QuiescentPeriodsResults{k});
                % Max
                QuiescentPeriodsStats{4,k+1} = max(QuiescentPeriodsResults{k});
                % Min
                QuiescentPeriodsStats{5,k+1} = min(QuiescentPeriodsResults{k});
                % N
                QuiescentPeriodsStats{6,k+1} = rowCounts(k);
            else
                rowCounts(k) = 0;
            end
        end

        QuiescentPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
        for k = 1:numel(sliceCells)
            QuiescentPeriodsOutput(1:rowCounts(k),k) = num2cell(QuiescentPeriodsResults{k});
        end
        TotalPeriods = reshape(QuiescentPeriodsOutput, [numel(QuiescentPeriodsOutput) ,1]);
        TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
        QuiescentPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));

        % Mean
        QuiescentPeriodsStats{1,end} = mean(TotalPeriods);
        % SE
        QuiescentPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
        % Median
        QuiescentPeriodsStats{3,end} = median(TotalPeriods);
        % Max
        QuiescentPeriodsStats{4,end} = max(TotalPeriods);
        % Min
        QuiescentPeriodsStats{5,end} = min(TotalPeriods);
        % N
        QuiescentPeriodsStats{6,end} = numel(TotalPeriods);

        QuiescentPeriodsOutput = [[{'Descriptive Stats'}, QuiescentPeriodsHeaders, {'Quiescent Periods Window1'}]; QuiescentPeriodsStats; [{''}, QuiescentPeriodsHeaders, {'Quiescent Periods Window1'}];[cell(numel(TotalPeriods),1), QuiescentPeriodsOutput]];
        if ispc
            xlswrite(sliceFile, QuiescentPeriodsOutput, 'Quiescent Periods Window1');
        else
            xlwrite(sliceFile, QuiescentPeriodsOutput, 'Quiescent Periods Window1');
        end
        
        % Window 2
        QuiescentPeriodsHeaders = cell(1,numel(sliceCells));
        QuiescentPeriodsResults = cell(1,numel(sliceCells));
        rowCounts = zeros(1,numel(sliceCells));
    
        % Calculate stats
        QuiescentPeriodsStats = cell(6, numel(sliceCells) + 2);
        QuiescentPeriodsStats{1,1} = 'mean';
        QuiescentPeriodsStats{2,1} = 'SE';
        QuiescentPeriodsStats{3,1} = 'median';
        QuiescentPeriodsStats{4,1} = 'max';
        QuiescentPeriodsStats{5,1} = 'min';
        QuiescentPeriodsStats{6,1} = 'N';

        for k = 1:numel(sliceCells)
            QuiescentPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            if ~isempty(EventData(sliceCells(k)).QuiescentPeriods_Windows)
                QuiescentPeriodsResults{k} = [EventData(sliceCells(k)).QuiescentPeriods_Windows.Window2];
                rowCounts(k)  = numel([EventData(sliceCells(k)).QuiescentPeriods_Windows.Window2]);

                % Mean
                QuiescentPeriodsStats{1,k+1} = mean(QuiescentPeriodsResults{k});
                % SE
                QuiescentPeriodsStats{2,k+1} = std(QuiescentPeriodsResults{k})./sqrt(rowCounts(k));
                % Median
                QuiescentPeriodsStats{3,k+1} = median(QuiescentPeriodsResults{k});
                % Max
                QuiescentPeriodsStats{4,k+1} = max(QuiescentPeriodsResults{k});
                % Min
                QuiescentPeriodsStats{5,k+1} = min(QuiescentPeriodsResults{k});
                % N
                QuiescentPeriodsStats{6,k+1} = rowCounts(k);
            else
                rowCounts(k) = 0;
            end
        end

        QuiescentPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
        for k = 1:numel(sliceCells)
            QuiescentPeriodsOutput(1:rowCounts(k),k) = num2cell(QuiescentPeriodsResults{k});
        end
        TotalPeriods = reshape(QuiescentPeriodsOutput, [numel(QuiescentPeriodsOutput) ,1]);
        TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
        QuiescentPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));

        % Mean
        QuiescentPeriodsStats{1,end} = mean(TotalPeriods);
        % SE
        QuiescentPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
        % Median
        QuiescentPeriodsStats{3,end} = median(TotalPeriods);
        % Max
        QuiescentPeriodsStats{4,end} = max(TotalPeriods);
        % Min
        QuiescentPeriodsStats{5,end} = min(TotalPeriods);
        % N
        QuiescentPeriodsStats{6,end} = numel(TotalPeriods);

        QuiescentPeriodsOutput = [[{'Descriptive Stats'}, QuiescentPeriodsHeaders, {'Quiescent Periods Window2'}]; QuiescentPeriodsStats; [{''}, QuiescentPeriodsHeaders, {'Quiescent Periods Window2'}];[cell(numel(TotalPeriods),1), QuiescentPeriodsOutput]];
        
        if ispc
            xlswrite(sliceFile, QuiescentPeriodsOutput, 'Quiescent Periods Window2');
        else
            xlwrite(sliceFile, QuiescentPeriodsOutput, 'Quiescent Periods Window2');
        end
    end










end