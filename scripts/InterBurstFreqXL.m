function InterBurstFreqXL(EventData, currSlice, sliceCells, sliceFile)

    InterBurstFreqHeaders = cell(1,numel(sliceCells));
    InterBurstFreqResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    InterBurstFreqStats = cell(5, numel(sliceCells) + 2);
    InterBurstFreqStats{1,1} = 'mean';
    InterBurstFreqStats{2,1} = 'SE';
    InterBurstFreqStats{3,1} = 'median';
    InterBurstFreqStats{4,1} = 'max';
    InterBurstFreqStats{5,1} = 'min';

    for k = 1:numel(sliceCells)
        InterBurstFreqHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).burstFreq)
            InterBurstFreqResults{k} = [EventData(sliceCells(k)).burstFreq.inter];
            rowCounts(k) = numel([EventData(sliceCells(k)).burstFreq.inter]);
            
            % Mean
            InterBurstFreqStats{1,k+1} = mean(InterBurstFreqResults{k});
            % SE
            InterBurstFreqStats{2,k+1} = std(InterBurstFreqResults{k})./sqrt(rowCounts(k));
            % Median
            InterBurstFreqStats{3,k+1} = median(InterBurstFreqResults{k});
            % Max
            InterBurstFreqStats{4,k+1} = max(InterBurstFreqResults{k});
            % Min
            InterBurstFreqStats{5,k+1} = min(InterBurstFreqResults{k});
        else
            rowCounts(k) = 0;
        end 
    end

    InterBurstFreqOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        InterBurstFreqOutput(1:rowCounts(k),k) = num2cell(InterBurstFreqResults{k});
    end
    TotalInterBurstFreq = reshape(InterBurstFreqOutput, [numel(InterBurstFreqOutput),1]);
    TotalInterBurstFreq = reshape(cell2mat(TotalInterBurstFreq), [sum(rowCounts),1]);
    InterBurstFreqOutput(:,end) = num2cell(sort(TotalInterBurstFreq));
    
    % Mean
    InterBurstFreqStats{1,end} = mean(TotalInterBurstFreq);
    % SE
    InterBurstFreqStats{2,end} = std(TotalInterBurstFreq)./sqrt(numel(TotalInterBurstFreq));
    % Median
    InterBurstFreqStats{3,end} = median(TotalInterBurstFreq);
    % Max
    InterBurstFreqStats{4,end} = max(TotalInterBurstFreq);
    % Min
    InterBurstFreqStats{5,end} = min(TotalInterBurstFreq);
    
    InterBurstFreqOutput = [[{'Descriptive Stats'}, InterBurstFreqHeaders, {'All Inter-Burst Frequencies'}]; InterBurstFreqStats; [{''}, InterBurstFreqHeaders, {'All Inter-Burst Frequencies'}]; [cell(numel(TotalInterBurstFreq),1), InterBurstFreqOutput]];

    if ispc
        xlswrite(sliceFile, InterBurstFreqOutput, 'InterBurstFreq');
    else
        xlwrite(sliceFile, InterBurstFreqOutput, 'InterBurstFreq');
    end





end