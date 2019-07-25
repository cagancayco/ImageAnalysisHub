function BurstTimingXL(EventData, currSlice, sliceCells, sliceFile)

BurstTimingHeaders = cell(1, numel(sliceCells));
BurstTimingValues = cell(1, numel(sliceCells));
rowCounts = zeros(1,numel(sliceCells));

for k = 1:numel(sliceCells)
    BurstTimingHeaders{k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
    if ~isempty(EventData(sliceCells(k)).burstTiming)
        BurstTimingValues{1,k} = [EventData(sliceCells(k)).burstTiming];
        rowCounts(k) = numel([EventData(sliceCells(k)).burstTiming]);
    else
        BurstTimingValues{1,k} = NaN;
        rowCounts(k) = 0;
    end
end

BurstTimingOutput = cell(max(rowCounts),numel(sliceCells));
for k = 1:numel(sliceCells)
    BurstTimingOutput(1:rowCounts(k),k) = num2cell(BurstTimingValues{k});
end

BurstTimingOutput = [BurstTimingHeaders; BurstTimingOutput];

if ispc
    xlswrite(sliceFile, BurstTimingOutput, 'BurstTiming');
else
    xlwrite(sliceFile, BurstTimingOutput, 'BurstTiming');
end

end