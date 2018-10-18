function BurstDurationXL(EventData, currSlice, sliceCells, sliceFile)

    BurstDurHeaders = cell(1,numel(sliceCells));
    BurstDurResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    BurstDurStats = cell(5, numel(sliceCells) + 2);
    BurstDurStats{1,1} = 'mean';
    BurstDurStats{2,1} = 'SE';
    BurstDurStats{3,1} = 'median';
    BurstDurStats{4,1} = 'max';
    BurstDurStats{5,1} = 'min';

    for k = 1:numel(sliceCells)
        BurstDurHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).burstDuration)
            BurstDurResults{k} = [EventData(sliceCells(k)).burstDuration.duration];
            rowCounts(k) = numel([EventData(sliceCells(k)).burstDuration.duration]);
            
            % Mean
            BurstDurStats{1,k+1} = mean(BurstDurResults{k});
            % SE
            BurstDurStats{2,k+1} = std(BurstDurResults{k})./sqrt(rowCounts(k));
            % Median
            BurstDurStats{3,k+1} = median(BurstDurResults{k});
            % Max
            BurstDurStats{4,k+1} = max(BurstDurResults{k});
            % Min
            BurstDurStats{5,k+1} = min(BurstDurResults{k});
        else
            rowCounts(k) = 0;
        end 
    end

    BurstDurOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        BurstDurOutput(1:rowCounts(k),k) = num2cell(BurstDurResults{k});
    end
    TotalBurstDur = reshape(BurstDurOutput, [numel(BurstDurOutput),1]);
    TotalBurstDur = reshape(cell2mat(TotalBurstDur), [sum(rowCounts),1]);
    BurstDurOutput(:,end) = num2cell(TotalBurstDur);
    
    % Mean
    BurstDurStats{1,end} = mean(TotalBurstDur);
    % SE
    BurstDurStats{2,end} = std(TotalBurstDur)./sqrt(numel(TotalBurstDur));
    % Median
    BurstDurStats{3,end} = median(TotalBurstDur);
    % Max
    BurstDurStats{4,end} = max(TotalBurstDur);
    % Min
    BurstDurStats{5,end} = min(TotalBurstDur);
    
    BurstDurOutput = [[{'Descriptive Stats'}, BurstDurHeaders, {'All Burst Durations'}]; BurstDurStats; [{''}, BurstDurHeaders, {'All Burst Durations'}]; [cell(numel(TotalBurstDur),1), BurstDurOutput]];

    if ispc
        xlswrite(sliceFile, BurstDurOutput, 'BurstDuration');
    else
        xlwrite(sliceFile, BurstDurOutput, 'BurstDuration');
    end
end