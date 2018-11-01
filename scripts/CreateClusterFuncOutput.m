function CreateClusterFuncOutput(EventData, currSlice, sliceCells, outputDir)

    data = cell(numel(sliceCells), 1);
   
    for i = 1:numel(data)
        data{i} = EventData(sliceCells(i)).EventTiming * 1000;
    end
   
    timeCmsec = [EventData(sliceCells(1)).time]'*1000;
    
    outputFolder = fullfile(outputDir,'clusterFunc');
    if ~exist(outputFolder); mkdir(outputFolder); end
    
    fileName = strcat(currSlice,'_clusterFunc.mat');
    fileName = fullfile(outputFolder, fileName);
    
    save(fileName, 'data', 'timeCmsec')
    


end