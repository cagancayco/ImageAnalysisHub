function runAllThresholding(app)
clc

now = string(datetime('now','Format','d_MMM_y_HH.mm'));
matFilename = strcat('PeriodThresholds_',now,'.mat');
matPath = fullfile(app.OutputDirectoryEditField.Value,'MAT-files',matFilename);
excelFilename = strcat('PeriodThresholds_',now,'.xlsx');
excelPath = fullfile(app.OutputDirectoryEditField.Value,excelFilename);

warning off
EventData = loadThresholdData(app);

nCells = length(EventData);
for i = 1:nCells
    currCell = EventData(i);

    try
        numEvents = length(EventTimes(currCell));
    catch
        numEvents = 0;
    end
    if numEvents < 3
        runCell = 0;
    else
        runCell = 1;
    end  

    switch runCell

        case 1
            EventData(i).IEI = AllPeriods(currCell);
            EventData(i).IEI = EventData(i).IEI(2:end-1);
        case 0     
    end
end

headers = {'Data Source', 'Truncated Gaussian Single', 'log Truncated Gaussian Single', ...
                          'Gaussian Only Intersection', 'log Gaussian Only Intersection', ...
                          'Gaussian Only Minimum', 'log Gaussian Only Minimum', ...
                          'Gaussian Exponential Intersection', 'log Gaussian Exp-Exponential Intersection', ...
                          'Gaussian Exponential Minimum', 'log Gaussian Exp-Exponential Minimum',...
                          'Bodova'};

    
%% Run Pooled Thresholding

periods = {EventData.IEI}';
periods = cellfun(@(X) X, periods, 'UniformOutput', 0);
periods = cell2mat(periods(~cellfun('isempty', periods)));

pooled_thresholds = cell(1,12);

pooled_thresholds{1,1} = 'All Data';
pooled_thresholds{1,2} = thTruncGaussianSingle(periods, 0);
pooled_thresholds{1,3} = thTruncGaussianSingle(log(periods), 1);
[pooled_thresholds{1,6},pooled_thresholds{1,4}] = thGaussianOnly(periods, 0);
[pooled_thresholds{1,7},pooled_thresholds{1,5}] = thGaussianOnly(log(periods), 1);
[pooled_thresholds{1,10},pooled_thresholds{1,8}] = thGaussExp(periods);
[pooled_thresholds{1,11},pooled_thresholds{1,9}] = thGaussExpExp(log(periods));
pooled_thresholds{1,12} = thBodova(periods);

%% Run Condition Thresholding

slices = {EventData.sliceName};
idx = num2cell(cellfun(@(X) strfind(X, '_data'), slices));
[conditions, ~, condIdx] = unique(cellfun(@(X,Y) X(1:Y), slices, idx, 'UniformOutput', 0));

conditionsPeriods = cell(numel(conditions),1);

for i = 1:numel(condIdx)
    conditionsPeriods{condIdx(i)} = [conditionsPeriods{condIdx(i)}; EventData(i).IEI];
end

condition_thresholds = cell(numel(conditions),12);

for j = 1:numel(conditions)
    condition_thresholds{j,1} = conditions{j};
    condition_thresholds{j,2} = thTruncGaussianSingle(conditionsPeriods{j}, 0);
    condition_thresholds{j,3} = thTruncGaussianSingle(log(conditionsPeriods{j}), 1);
    [condition_thresholds{j,6},condition_thresholds{j,4}] = thGaussianOnly(conditionsPeriods{j}, 0);
    [condition_thresholds{j,7},condition_thresholds{j,5}] = thGaussianOnly(log(conditionsPeriods{j}), 1);
    [condition_thresholds{j,10},condition_thresholds{j,8}] = thGaussExp(conditionsPeriods{j});
    [condition_thresholds{j,11},condition_thresholds{j,9}] = thGaussExpExp(log(conditionsPeriods{j}));
    condition_thresholds{j,12} = thBodova(conditionsPeriods{j});    
end

%% Run Experiment Thresholding

[slices, ~, sliceIdx] = unique({EventData.sliceName});
    
slicesPeriods = cell(numel(slices),1);
for i = 1:numel(sliceIdx)
    slicesPeriods{sliceIdx(i)} = [slicesPeriods{sliceIdx(i)}; EventData(i).IEI];
end

slice_thresholds = cell(numel(slices),12);

for j = 1:numel(slices)
    slice_thresholds{j,1} = slices{j};
    slice_thresholds{j,2} = thTruncGaussianSingle(slicesPeriods{j}, 0);
    slice_thresholds{j,3} = thTruncGaussianSingle(log(slicesPeriods{j}), 1);
    [slice_thresholds{j,6},slice_thresholds{j,4}] = thGaussianOnly(slicesPeriods{j}, 0);
    [slice_thresholds{j,7},slice_thresholds{j,5}] = thGaussianOnly(log(slicesPeriods{j}), 1);
    [slice_thresholds{j,10},slice_thresholds{j,8}] = thGaussExp(slicesPeriods{j});
    [slice_thresholds{j,11},slice_thresholds{j,9}] = thGaussExpExp(log(slicesPeriods{j}));
    slice_thresholds{j,12} = thBodova(slicesPeriods{j});  
end


%% Run Cell Thresholding

cell_thresholds = cell(length(EventData),12);

for j = 1:length(EventData)
    cell_thresholds{j,1} = EventData(j).cellName;
    if ~isempty(EventData(j).IEI)
        cell_thresholds{j,2} = thTruncGaussianSingle(EventData(j).IEI, 0);
        cell_thresholds{j,3} = thTruncGaussianSingle(log(EventData(j).IEI), 1);
        [cell_thresholds{j,6},cell_thresholds{j,4}] = thGaussianOnly(EventData(j).IEI, 0);
        [cell_thresholds{j,7},cell_thresholds{j,5}] = thGaussianOnly(log(EventData(j).IEI), 1);
        [cell_thresholds{j,10},cell_thresholds{j,8}] = thGaussExp(EventData(j).IEI);
        [cell_thresholds{j,11},cell_thresholds{j,9}] = thGaussExpExp(log(EventData(j).IEI));
        cell_thresholds{j,12} = thBodova(EventData(j).IEI);  
        
    end
end



%% Combine results and write to Excel and MAT-file


thresholds = [headers; pooled_thresholds; condition_thresholds; slice_thresholds; cell_thresholds];
save(matPath, 'thresholds')

if ispc
    xlswrite(excelPath, thresholds, 'Thresholds');
else    
    xlwrite(excelPath, thresholds, 'Thresholds');
end


end