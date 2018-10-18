function InterBurstIntXL(EventData, currSlice, sliceCells, sliceFile)

    InterBurstIntHeaders = cell(1,numel(sliceCells));
    InterBurstIntResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    InterBurstIntStats = cell(5, numel(sliceCells) + 2);
    InterBurstIntStats{1,1} = 'mean';
    InterBurstIntStats{2,1} = 'SE';
    InterBurstIntStats{3,1} = 'median';
    InterBurstIntStats{4,1} = 'max';
    InterBurstIntStats{5,1} = 'min';

    for k = 1:numel(sliceCells)
        InterBurstIntHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).interburstInt)
            InterBurstIntResults{k} = [EventData(sliceCells(k)).interburstInt];
            rowCounts(k) = numel([EventData(sliceCells(k)).interburstInt]);
            
            % Mean
            InterBurstIntStats{1,k+1} = mean(InterBurstIntResults{k});
            % SE
            InterBurstIntStats{2,k+1} = std(InterBurstIntResults{k})./sqrt(rowCounts(k));
            % Median
            InterBurstIntStats{3,k+1} = median(InterBurstIntResults{k});
            % Max
            InterBurstIntStats{4,k+1} = max(InterBurstIntResults{k});
            % Min
            InterBurstIntStats{5,k+1} = min(InterBurstIntResults{k});
        else
            rowCounts(k) = 0;
        end 
    end

    InterBurstIntOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        InterBurstIntOutput(1:rowCounts(k),k) = num2cell(InterBurstIntResults{k});
    end
    TotalInterBurstInt = reshape(InterBurstIntOutput, [numel(InterBurstIntOutput),1]);
    TotalInterBurstInt = reshape(cell2mat(TotalInterBurstInt), [sum(rowCounts),1]);
    InterBurstIntOutput(:,end) = num2cell(sort(TotalInterBurstInt));
    
    % Mean
    InterBurstIntStats{1,end} = mean(TotalInterBurstInt);
    % SE
    InterBurstIntStats{2,end} = std(TotalInterBurstInt)./sqrt(numel(TotalInterBurstInt));
    % Median
    InterBurstIntStats{3,end} = median(TotalInterBurstInt);
    % Max
    InterBurstIntStats{4,end} = max(TotalInterBurstInt);
    % Min
    InterBurstIntStats{5,end} = min(TotalInterBurstInt);
    
    InterBurstIntOutput = [[{'Descriptive Stats'}, InterBurstIntHeaders, {'All Inter-Burst Intervals'}]; InterBurstIntStats; [{''}, InterBurstIntHeaders, {'All Inter-Burst Intervals'}]; [cell(numel(TotalInterBurstInt),1), InterBurstIntOutput]];

    if ispc
        xlswrite(sliceFile, InterBurstIntOutput, 'InterBurstInt');
    else
        xlwrite(sliceFile, InterBurstIntOutput, 'InterBurstInt');
    end





end