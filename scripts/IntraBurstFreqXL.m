function IntraBurstFreqXL(EventData, currSlice, sliceCells, sliceFile)

    IntraBurstFreqHeaders = cell(1,numel(sliceCells));
    IntraBurstFreqResults = cell(1,numel(sliceCells));
    rowCounts = zeros(1,numel(sliceCells));
    
    % Calculate stats
    IntraBurstFreqStats = cell(5, numel(sliceCells) + 2);
    IntraBurstFreqStats{1,1} = 'mean';
    IntraBurstFreqStats{2,1} = 'SE';
    IntraBurstFreqStats{3,1} = 'median';
    IntraBurstFreqStats{4,1} = 'max';
    IntraBurstFreqStats{5,1} = 'min';

    for k = 1:numel(sliceCells)
        IntraBurstFreqHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
        if ~isempty(EventData(sliceCells(k)).burstFreq)
            IntraBurstFreqResults{k} = [EventData(sliceCells(k)).burstFreq.intra.frequency];
            rowCounts(k) = numel([EventData(sliceCells(k)).burstFreq.intra.frequency]);
            
            % Mean
            IntraBurstFreqStats{1,k+1} = mean(IntraBurstFreqResults{k});
            % SE
            IntraBurstFreqStats{2,k+1} = std(IntraBurstFreqResults{k})./sqrt(rowCounts(k));
            % Median
            IntraBurstFreqStats{3,k+1} = median(IntraBurstFreqResults{k});
            % Max
            IntraBurstFreqStats{4,k+1} = max(IntraBurstFreqResults{k});
            % Min
            IntraBurstFreqStats{5,k+1} = min(IntraBurstFreqResults{k});
        else
            rowCounts(k) = 0;
        end 
    end

    IntraBurstFreqOutput = cell(sum(rowCounts),numel(sliceCells)+1);
    for k = 1:numel(sliceCells)
        IntraBurstFreqOutput(1:rowCounts(k),k) = num2cell(IntraBurstFreqResults{k});
    end
    TotalIntraBurstFreq = reshape(IntraBurstFreqOutput, [numel(IntraBurstFreqOutput),1]);
    TotalIntraBurstFreq = reshape(cell2mat(TotalIntraBurstFreq), [sum(rowCounts),1]);
    IntraBurstFreqOutput(:,end) = num2cell(sort(TotalIntraBurstFreq));
    
    % Mean
    IntraBurstFreqStats{1,end} = mean(TotalIntraBurstFreq);
    % SE
    IntraBurstFreqStats{2,end} = std(TotalIntraBurstFreq)./sqrt(numel(TotalIntraBurstFreq));
    % Median
    IntraBurstFreqStats{3,end} = median(TotalIntraBurstFreq);
    % Max
    IntraBurstFreqStats{4,end} = max(TotalIntraBurstFreq);
    % Min
    IntraBurstFreqStats{5,end} = min(TotalIntraBurstFreq);
    
    IntraBurstFreqOutput = [[{'Descriptive Stats'}, IntraBurstFreqHeaders, {'All Intra-Burst Frequencies'}]; IntraBurstFreqStats; [{''}, IntraBurstFreqHeaders, {'All Intra-Burst Frequencies'}]; [cell(numel(TotalIntraBurstFreq),1), IntraBurstFreqOutput]];

    if ispc
        xlswrite(sliceFile, IntraBurstFreqOutput, 'IntraBurstFreq');
    else
        xlwrite(sliceFile, IntraBurstFreqOutput, 'IntraBurstFreq');
    end
    







end