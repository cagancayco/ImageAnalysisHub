function AllPeriodsXL(EventData, currSlice, sliceCells, sliceFile)
    AllPeriodsHeaders = cell(1,numel(sliceCells));
    AllPeriodsResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    AllPeriodsStats = cell(6, numel(sliceCells) + 2);
    AllPeriodsStats{1,1} = 'mean';
    AllPeriodsStats{2,1} = 'SE';
    AllPeriodsStats{3,1} = 'median';
    AllPeriodsStats{4,1} = 'max';
    AllPeriodsStats{5,1} = 'min';
    AllPeriodsStats{6,1} = 'N';
            
    for k = 1:numel(sliceCells)
        AllPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).AllPeriods)
            AllPeriodsResults{k} = [EventData(sliceCells(k)).AllPeriods];
            rowCounts(k)  = numel([EventData(sliceCells(k)).AllPeriods]);
            
            % Mean
            AllPeriodsStats{1,k+1} = mean(AllPeriodsResults{k});
            % SE
            AllPeriodsStats{2,k+1} = std(AllPeriodsResults{k})./sqrt(rowCounts(k));
            % Median
            AllPeriodsStats{3,k+1} = median(AllPeriodsResults{k});
            % Max
            AllPeriodsStats{4,k+1} = max(AllPeriodsResults{k});
            % Min
            AllPeriodsStats{5,k+1} = min(AllPeriodsResults{k});
            % N
            AllPeriodsStats{6,k+1} = rowCounts(k);
        else
            rowCounts(k) = 0;
        end
    end
            
    AllPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        AllPeriodsOutput(1:rowCounts(k),k) = num2cell(AllPeriodsResults{k});
    end
    TotalPeriods = reshape(AllPeriodsOutput, [numel(AllPeriodsOutput) ,1]);
    TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
    AllPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));
%     AllPeriodsList = num2cell(sort(TotalPeriods));
    
    % Mean
    AllPeriodsStats{1,end} = mean(TotalPeriods);
    % SE
    AllPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
    % Median
    AllPeriodsStats{3,end} = median(TotalPeriods);
    % Max
    AllPeriodsStats{4,end} = max(TotalPeriods);
    % Min
    AllPeriodsStats{5,end} = min(TotalPeriods);
    % N
    AllPeriodsStats{6,end} = numel(TotalPeriods);
            
    AllPeriodsOutput = [[{'Descriptive Stats'}, AllPeriodsHeaders, {'All Periods'}]; AllPeriodsStats; [{''}, AllPeriodsHeaders, {'All Periods'}];[cell(numel(TotalPeriods),1), AllPeriodsOutput]];
%     lenAllPeriodsOutput = size(AllPeriodsOutput,1);
%     AllPeriodsList = [AllPeriodsList; cell(lenAllPeriodsOutput - length(AllPeriodsList),1)];
%     AllPeriodsOutput = [AllPeriodsOutput, cell(lenAllPeriodsOutput,1), AllPeriodsList];
    
    if ispc
        xlswrite(sliceFile, AllPeriodsOutput, 'All Periods');
    else
        xlwrite(sliceFile, AllPeriodsOutput, 'All Periods');
    end
            
%% Window Analysis

    if ~isempty(EventData(1).analysisWindows)
        % Window 1
        AllPeriodsHeaders = cell(1,numel(sliceCells));
        AllPeriodsResults = cell(1,numel(sliceCells));
        rowCounts = zeros(1,numel(sliceCells));
    
        % Calculate stats
        AllPeriodsStats = cell(6, numel(sliceCells) + 2);
        AllPeriodsStats{1,1} = 'mean';
        AllPeriodsStats{2,1} = 'SE';
        AllPeriodsStats{3,1} = 'median';
        AllPeriodsStats{4,1} = 'max';
        AllPeriodsStats{5,1} = 'min';
        AllPeriodsStats{6,1} = 'N';
        

        for k = 1:numel(sliceCells)
            AllPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            if ~isempty(EventData(sliceCells(k)).AllPeriods_Windows)
                AllPeriodsResults{k} = [EventData(sliceCells(k)).AllPeriods_Windows.Window1];
                rowCounts(k)  = numel([EventData(sliceCells(k)).AllPeriods_Windows.Window1]);

                % Mean
                AllPeriodsStats{1,k+1} = mean(AllPeriodsResults{k});
                % SE
                AllPeriodsStats{2,k+1} = std(AllPeriodsResults{k})./sqrt(rowCounts(k));
                % Median
                AllPeriodsStats{3,k+1} = median(AllPeriodsResults{k});
                % Max
                AllPeriodsStats{4,k+1} = max(AllPeriodsResults{k});
                % Min
                AllPeriodsStats{5,k+1} = min(AllPeriodsResults{k});
                % N
                AllPeriodsStats{6,k+1} = rowCounts(k);
            else
                rowCounts(k) = 0;
            end
        end

        AllPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
        for k = 1:numel(sliceCells)
            AllPeriodsOutput(1:rowCounts(k),k) = num2cell(AllPeriodsResults{k});
        end
        TotalPeriods = reshape(AllPeriodsOutput, [numel(AllPeriodsOutput) ,1]);
        TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
        AllPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));
%         AllPeriodsList = num2cell(sort(TotalPeriods));
        % Mean
        AllPeriodsStats{1,end} = mean(TotalPeriods);
        % SE
        AllPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
        % Median
        AllPeriodsStats{3,end} = median(TotalPeriods);
        % Max
        AllPeriodsStats{4,end} = max(TotalPeriods);
        % Min
        AllPeriodsStats{5,end} = min(TotalPeriods);
        % N
        AllPeriodsStats{6,end} = numel(TotalPeriods);

        AllPeriodsOutput = [[{'Descriptive Stats'}, AllPeriodsHeaders, {'All Periods Window1'}]; AllPeriodsStats; [{''}, AllPeriodsHeaders, {'All Periods Window1'}];[cell(numel(TotalPeriods),1), AllPeriodsOutput]];
%         lenAllPeriodsOutput = size(AllPeriodsOutput,1);
%         AllPeriodsList = [AllPeriodsList; cell(lenAllPeriodsOutput - length(AllPeriodsList),1)];
%         AllPeriodsOutput = [AllPeriodsOutput, cell(lenAllPeriodsOutput,1), AllPeriodsList];
        if ispc
            xlswrite(sliceFile, AllPeriodsOutput, 'All Periods Window1');
        else
            xlwrite(sliceFile, AllPeriodsOutput, 'All Periods Window1');
        end
        
        % Window 2
        AllPeriodsHeaders = cell(1,numel(sliceCells));
        AllPeriodsResults = cell(1,numel(sliceCells));
        rowCounts = zeros(1,numel(sliceCells));
    
        % Calculate stats
        AllPeriodsStats = cell(6, numel(sliceCells) + 2);
        AllPeriodsStats{1,1} = 'mean';
        AllPeriodsStats{2,1} = 'SE';
        AllPeriodsStats{3,1} = 'median';
        AllPeriodsStats{4,1} = 'max';
        AllPeriodsStats{5,1} = 'min';
        AllPeriodsStats{6,1} = 'N';

        for k = 1:numel(sliceCells)
            AllPeriodsHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            if ~isempty(EventData(sliceCells(k)).AllPeriods_Windows)
                AllPeriodsResults{k} = [EventData(sliceCells(k)).AllPeriods_Windows.Window2];
                rowCounts(k)  = numel([EventData(sliceCells(k)).AllPeriods_Windows.Window2]);

                % Mean
                AllPeriodsStats{1,k+1} = mean(AllPeriodsResults{k});
                % SE
                AllPeriodsStats{2,k+1} = std(AllPeriodsResults{k})./sqrt(rowCounts(k));
                % Median
                AllPeriodsStats{3,k+1} = median(AllPeriodsResults{k});
                % Max
                AllPeriodsStats{4,k+1} = max(AllPeriodsResults{k});
                % Min
                AllPeriodsStats{5,k+1} = min(AllPeriodsResults{k});
                % N
                AllPeriodsStats{6,k+1} = rowCounts(k);
            else
                rowCounts(k) = 0;
            end
        end

        AllPeriodsOutput = cell(sum(rowCounts),numel(sliceCells)+1);
        for k = 1:numel(sliceCells)
            AllPeriodsOutput(1:rowCounts(k),k) = num2cell(AllPeriodsResults{k});
        end
        TotalPeriods = reshape(AllPeriodsOutput, [numel(AllPeriodsOutput) ,1]);
        TotalPeriods = reshape(cell2mat(TotalPeriods), [sum(rowCounts), 1]);
        AllPeriodsOutput(:,end) = num2cell(sort(TotalPeriods));
%         AllPeriodsList = num2cell(sort(TotalPeriods));
        % Mean
        AllPeriodsStats{1,end} = mean(TotalPeriods);
        % SE
        AllPeriodsStats{2,end} = std(TotalPeriods)./sqrt(numel(TotalPeriods));
        % Median
        AllPeriodsStats{3,end} = median(TotalPeriods);
        % Max
        AllPeriodsStats{4,end} = max(TotalPeriods);
        % Min
        AllPeriodsStats{5,end} = min(TotalPeriods);
        % N
        AllPeriodsStats{6,end} = numel(TotalPeriods);

        AllPeriodsOutput = [[{'Descriptive Stats'}, AllPeriodsHeaders, {'All Periods Window2'}]; AllPeriodsStats; [{''}, AllPeriodsHeaders, {'All Periods Window2'}];[cell(numel(TotalPeriods),1), AllPeriodsOutput]];
%         lenAllPeriodsOutput = size(AllPeriodsOutput,1);
%         AllPeriodsList = [AllPeriodsList; cell(lenAllPeriodsOutput - length(AllPeriodsList),1)];
%         AllPeriodsOutput = [AllPeriodsOutput, cell(lenAllPeriodsOutput,1), AllPeriodsList];
        if ispc
            xlswrite(sliceFile, AllPeriodsOutput, 'All Periods Window2');
        else
            xlwrite(sliceFile, AllPeriodsOutput, 'All Periods Window2');
        end
    end






end