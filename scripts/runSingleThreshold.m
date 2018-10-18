function EventData = runSingleThreshold(app, EventData)

plotThresholdFit = app.PlotfitwithPreviewCheckBox.Value;

if app.PooledButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
    thresholdGroup = 'Pooled_';
    periods = {EventData.IEI}';
    periods = cellfun(@(X) X, periods, 'UniformOutput', 0);
    periods = cell2mat(periods(~cellfun('isempty', periods)));
    
    if app.TruncatedGaussianSingleButton.Value
        threshold = thTruncGaussianSingle(periods, 0);
        thresholdMethod = [thresholdGroup,'TruncatedGaussianSingle'];
    elseif app.logTruncatedGaussianSingleButton.Value
        threshold = thTruncGaussianSingle(log(periods), 1);
        thresholdMethod = [thresholdGroup,'logTrucatedGaussianSingle'];
    elseif app.GaussianOnlyIntersectionButton.Value
        [~, threshold] = thGaussianOnly(periods, 0);
        thresholdMethod = [thresholdGroup,'GaussianOnlyIntersection'];
    elseif app.logGaussianOnlyIntersectionButton.Value
        [~, threshold] = thGaussianOnly(log(periods), 1);
        thresholdMethod = [thresholdGroup,'logGaussianOnlyIntersection'];
    elseif app.GaussianOnlyMinimumButton.Value
        [threshold, ~] = thGaussianOnly(periods, 0);
        thresholdMethod = [thresholdGroup,'GaussianOnlyMinimum'];
    elseif app.logGaussianOnlyMinimumButton.Value
        [threshold, ~] = thGaussianOnly(log(periods), 1);  
        thresholdMethod = [thresholdGroup,'logGaussianOnlyMinimum'];
    elseif app.GaussianExponentialIntersectionButton.Value
        [~, threshold] = thGaussExp(periods);
        thresholdMethod = [thresholdGroup,'GaussianExponentialIntersection'];
    elseif app.logGaussianExpExponentialIntersectionButton.Value
        [~, threshold] = thGaussExpExp(log(periods));
        thresholdMethod = [thresholdGroup,'logGaussianExpExponentialIntersection'];
    elseif app.GaussianExponentialMinimumButton.Value
        [threshold, ~] = thGaussExp(periods);
        thresholdMethod = [thresholdGroup,'GaussianExponentialMinimum'];
    elseif app.logGaussianExpExponentialMinimumButton.Value
        [threshold, ~] = thGaussExpExp(log(periods));    
        thresholdMethod = [thresholdGroup,'logGaussianExpExponentialMinimum'];
    elseif app.BodovaButton.Value
        threshold = thBodova(periods);
        thresholdMethod = [thresholdGroup,'Bodova'];
    end
    
    for i = 1:length(EventData)
        EventData(i).thresholdMethod = thresholdMethod;
        EventData(i).threshold = threshold;
    end
     
        
    
    
elseif app.ConditionLevelButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
    thresholdGroup = 'Condition_';
    slices = {EventData.sliceName};
    idx = num2cell(cellfun(@(X) strfind(X, '_data'), slices));
    [conditions, ~, condIdx] = unique(cellfun(@(X,Y) X(1:Y), slices, idx, 'UniformOutput', 0));
    
    conditionsPeriods = cell(numel(conditions),1);
    
    for i = 1:numel(condIdx)
        conditionsPeriods{condIdx(i)} = [conditionsPeriods{condIdx(i)}; EventData(i).IEI];
    end
    
    
    threshold = zeros(numel(conditions),1);
    
    if app.TruncatedGaussianSingleButton.Value
        for j = 1:numel(conditions)
            threshold(j) = thTruncGaussianSingle(conditionsPeriods{j}, 0);
        end
        thresholdMethod = [thresholdGroup,'TruncatedGaussianSingle'];
    elseif app.logTruncatedGaussianSingleButton.Value
        for j = 1:numel(conditions)
            threshold(j) = thTruncGaussianSingle(log(conditionsPeriods{j}), 1);
        end
        thresholdMethod = [thresholdGroup,'logTrucatedGaussianSingle'];
    elseif app.GaussianOnlyIntersectionButton.Value
        for j = 1:numel(conditions)
            [~, threshold(j)] = thGaussianOnly(conditionsPeriods{j}, 0);
        end
        thresholdMethod = [thresholdGroup,'GaussianOnlyIntersection'];
    elseif app.logGaussianOnlyIntersectionButton.Value
        for j = 1:numel(conditions)
            [~, threshold(j)] = thGaussianOnly(log(conditionsPeriods{j}), 1);
        end
        thresholdMethod = [thresholdGroup,'logGaussianOnlyIntersection'];
    elseif app.GaussianOnlyMinimumButton.Value
        for j = 1:numel(conditions)
            [threshold(j), ~] = thGaussianOnly(conditionsPeriods{j}, 0);
        end
        thresholdMethod = [thresholdGroup,'GaussianOnlyMinimum'];
    elseif app.logGaussianOnlyMinimumButton.Value
        for j = 1:numel(conditions)
            [threshold(j), ~] = thGaussianOnly(log(conditionsPeriods{j}), 1);
        end 
        thresholdMethod = [thresholdGroup,'logGaussianOnlyMinimum'];
    elseif app.GaussianExponentialIntersectionButton.Value
        for j = 1:numel(conditions)
            [~, threshold(j)] = thGaussExp(conditionsPeriods{j});
        end
        thresholdMethod = [thresholdGroup,'GaussianExponentialIntersection'];
    elseif app.logGaussianExpExponentialIntersectionButton.Value
        for j = 1:numel(conditions)
            [~, threshold(j)] = thGaussExpExp(log(conditionsPeriods{j}));
        end
        thresholdMethod = [thresholdGroup,'logGaussuanExpExponentialIntersection'];
    elseif app.GaussianExponentialMinimumButton.Value
        for j = 1:numel(conditions)
            [threshold(j), ~] = thGaussExp(conditionsPeriods{j});
        end
        thresholdMethod = [thresholdGroup,'GaussianExponentialMinimum'];
    elseif app.logGaussianExpExponentialMinimumButton.Value
        for j = 1:numel(conditions)
            [threshold(j), ~] = thGaussExpExp(log(conditionsPeriods{j}));
        end
        thresholdMethod = [thresholdGroup,'logGaussianExpExponentialMinimum'];
    elseif app.BodovaButton.Value
        for j = 1:numel(conditions)
            threshold(j) = thBodova(conditionsPeriods{j});
        end
        thresholdMethod = [thresholdGroup,'Bodova'];
    end
    
    
    for k = 1:length(EventData)
        EventData(k).thresholdMethod = thresholdMethod;
        EventData(k).threshold = threshold(condIdx(k));
    end
    
elseif app.ExperimentLevelButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
    thresholdGroup = 'Experiment_';
    [slices, ~, sliceIdx] = unique({EventData.sliceName});
    
    slicesPeriods = cell(numel(slices),1);
    for i = 1:numel(sliceIdx)
        slicesPeriods{sliceIdx(i)} = [slicesPeriods{sliceIdx(i)}; EventData(i).IEI];
    end
    
    threshold = zeros(numel(slices),1);
    
    if app.TruncatedGaussianSingleButton.Value
        for j = 1:numel(slices)
            threshold(j) = thTruncGaussianSingle(slicesPeriods{j}, 0);
        end
        thresholdMethod = [thresholdGroup,'TruncatedGaussianSingle'];
    elseif app.logTruncatedGaussianSingleButton.Value
        for j = 1:numel(slices)
            threshold(j) = thTruncGaussianSingle(log(slicesPeriods{j}), 1);
        end
        thresholdMethod = [thresholdGroup,'logTrucatedGaussianSingle'];
    elseif app.GaussianOnlyIntersectionButton.Value
        for j = 1:numel(slices)
            [~, threshold(j)] = thGaussianOnly(slicesPeriods{j}, 0);
        end
        thresholdMethod = [thresholdGroup,'GaussianOnlyIntersection'];
    elseif app.logGaussianOnlyIntersectionButton.Value
        for j = 1:numel(cslices)
            [~, threshold(j)] = thGaussianOnly(log(slicesPeriods{j}), 1);
        end
        thresholdMethod = [thresholdGroup,'logGaussianOnlyIntersection'];
    elseif app.GaussianOnlyMinimumButton.Value
        for j = 1:numel(slices)
            [threshold(j), ~] = thGaussianOnly(slicesPeriods{j}, 0);
        end
        thresholdMethod = [thresholdGroup,'GaussianOnlyMinimum'];
    elseif app.logGaussianOnlyMinimumButton.Value
        for j = 1:numel(slices)
            [threshold(j), ~] = thGaussianOnly(log(slicesPeriods{j}), 1);
        end 
        thresholdMethod = [thresholdGroup,'logGaussianOnlyMinimum'];
    elseif app.GaussianExponentialIntersectionButton.Value
        for j = 1:numel(slices)
            [~, threshold(j)] = thGaussExp(slicesPeriods{j});
        end
        thresholdMethod = [thresholdGroup,'GaussianExponentialIntersection'];
    elseif app.logGaussianExpExponentialIntersectionButton.Value
        for j = 1:numel(slices)
            [~, threshold(j)] = thGaussExpExp(log(slicesPeriods{j}));
        end
        thresholdMethod = [thresholdGroup,'logGaussuanExpExponentialIntersection'];
    elseif app.GaussianExponentialMinimumButton.Value
        for j = 1:numel(slices)
            [threshold(j), ~] = thGaussExp(slicesPeriods{j});
        end
        thresholdMethod = [thresholdGroup,'GaussianExponentialMinimum'];
    elseif app.logGaussianExpExponentialMinimumButton.Value
        for j = 1:numel(slices)
            [threshold(j), ~] = thGaussExpExp(log(slicesPeriods{j}));
        end
        thresholdMethod = [thresholdGroup,'logGaussianExpExponentialMinimum'];
    elseif app.BodovaButton.Value
        for j = 1:numel(slices)
            threshold(j) = thBodova(slicesPeriods{j});
        end
        thresholdMethod = [thresholdGroup,'Bodova'];
    end
    
    
    for k = 1:length(EventData)
        EventData(k).thresholdMethod = thresholdMethod;
        EventData(k).threshold = threshold(sliceIdx(k));
    end
    
elseif app.CellLevelButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
    thresholdGroup = 'Cell_';
    
    for i = 1:length(EventData)
        
        if ~isempty(EventData(i).IEI)
            if app.TruncatedGaussianSingleButton.Value
                EventData(i).threshold = thTruncGaussianSingle(EventData(i).IEI, 0);
                EventData(i).thresholdMethod = [thresholdGroup,'TruncatedGaussianSingle'];
            elseif app.logTruncatedGaussianSingleButton.Value
                EventData(i).threshold = thTruncGaussianSingle(log(EventData(i).IEI), 1);
                EventData(i).thresholdMethod = [thresholdGroup,'logTrucatedGaussianSingle'];
            elseif app.GaussianOnlyIntersectionButton.Value
                [~, EventData(i).threshold] = thGaussianOnly(EventData(i).IEI, 0);
                EventData(i).thresholdMethod = [thresholdGroup,'GaussianOnlyIntersection'];
            elseif app.logGaussianOnlyIntersectionButton.Value
                [~, EventData(i).threshold] = thGaussianOnly(log(EventData(i).IEI), 1);
                EventData(i).thresholdMethod = [thresholdGroup,'logGaussianOnlyIntersection'];
            elseif app.GaussianOnlyMinimumButton.Value
                [EventData(i).threshold, ~] = thGaussianOnly(EventData(i).IEI, 0);
                EventData(i).thresholdMethod = [thresholdGroup,'GaussianOnlyMinimum'];
            elseif app.logGaussianOnlyMinimumButton.Value
                [EventData(i).threshold, ~] = thGaussianOnly(log(EventData(i).IEI), 1);  
                EventData(i).thresholdMethod = [thresholdGroup,'logGaussianOnlyMinimum'];
            elseif app.GaussianExponentialIntersectionButton.Value
                [~, EventData(i).threshold] = thGaussExp(EventData(i).IEI);
                EventData(i).thresholdMethod = [thresholdGroup,'GaussianExponentialIntersection'];
            elseif app.logGaussianExpExponentialIntersectionButton.Value
                [~, EventData(i).threshold] = thGaussExpExp(log(EventData(i).IEI));
                EventData(i).thresholdMethod = [thresholdGroup,'logGaussuanExpExponentialIntersection'];
            elseif app.GaussianExponentialMinimumButton.Value
                [EventData(i).threshold, ~] = thGaussExp(EventData(i).IEI);
                EventData(i).thresholdMethod = [thresholdGroup,'GaussianExponentialMinimum'];
            elseif app.logGaussianExpExponentialMinimumButton.Value
                [EventData(i).threshold, ~] = thGaussExpExp(log(EventData(i).IEI));    
                EventData(i).thresholdMethod = [thresholdGroup,'logGaussianExpExponentialMinimum'];
            elseif app.BodovaButton.Value
                EventData(i).threshold = thBodova(EventData(i).IEI);
                EventData(i).thresholdMethod = [thresholdGroup,'Bodova'];
            end
        
        
        end
    end

    
end







end