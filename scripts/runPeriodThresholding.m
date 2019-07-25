function runPeriodThresholding(app)

warning off
clc

if ~ispc
    xlwrite_java_path = [pwd, '/poi_library'];
    javaaddpath(xlwrite_java_path)
end

if isfolder('2_Period_Thresholding/Plots')
    outputDir = '2_Period_Thresholding/Plots/';
else
    outputDir = [pwd,'/'];
end

    
% Create Output filename timestamp appendix
now = string(datetime('now','Format','d_MMM_y_HH.mm'));
matFilename = strcat(app.LoadEventDataMATfileEditField.Value(1:end-4),'_',now,'.mat');
excelFilename = strcat('_PeriodThresholdingz_',now,'.xlsx');

% Load Event Data MAT-file
EventData = load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));
nCells = length(EventData.EventData);

newEventData = EventData.EventData;
slices = unique({newEventData.sliceName});

plotTh = app.ExportPlotsCheckBox.Value;

% Run AllPeriods

for i = 1:nCells
    warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
    currCell = newEventData(i);

    if length(currCell.events) < 3
        runCell = 0;
    else
        runCell = 1;
    end

    switch runCell

        case 1
            newEventData(i).IEI = AllPeriods(currCell);
            newEventData(i).IEI = newEventData(i).IEI(2:end-1);
        case 0     
    end
end


if app.PooledCheckBox.Value
    alldata = struct;
    thresholdGroup = 'Pooled_';
    periods = {EventData.IEI}';
    periods = cellfun(@(X) X, periods, 'UniformOutput', 0);
    periods = cell2mat(periods(~cellfun('isempty', periods)));
    
    if app.GaussianOnlyIntersectionCheckBox.Value
        [~, alldata.th_GaussianOnlyInt, ~, ~, ~, plot_GaussianOnlyInt] = thGaussianOnly2(periods, 0, plotTh);
        if ~isnan(plot_GaussianOnlyInt)
            ax = plot_GaussianOnlyInt.CurrentAxes;
            set(get(ax, 'title'), 'string', strcat(thresholdGroup,'GaussianOnlyIntersection'));
            saveas(plot_GaussianOnlyInt, strcat(outputDir, thresholdGroup,'GaussianOnlyIntersection_', now, '.png'))
        end
    end
    
    if app.logGaussianOnlyIntersectionCheckBox.Value
        [~, app.th_logGaussianOnlyInt, ~, ~, ~] = thGaussianOnly2(log(periods), 1, plotTh);
    end
    
    
end

if app.ConditionLevelButton.Value
    thresholdGroup = 'Condition_';
    slices = {EventData.sliceName};
    idx = num2cell(cellfun(@(X) strfind(X, '_data'), slices));
    [conditions, ~, condIdx] = unique(cellfun(@(X,Y) X(1:Y), slices, idx, 'UniformOutput', 0));
    
    conditionsPeriods = cell(numel(conditions),1);
    
    for i = 1:numel(condIdx)
        conditionsPeriods{condIdx(i)} = [conditionsPeriods{condIdx(i)}; EventData(i).IEI];
    end
end

if app.ExperimentLevelButton.Value
    thresholdGroup = 'Experiment_';
    [slices, ~, sliceIdx] = unique({EventData.sliceName});
    
    slicesPeriods = cell(numel(slices),1);
    for i = 1:numel(sliceIdx)
        slicesPeriods{sliceIdx(i)} = [slicesPeriods{sliceIdx(i)}; EventData(i).IEI];
    end
end

if app.CellLevelButton.Value
    thresholdGroup = 'Cell_';
    
end



end