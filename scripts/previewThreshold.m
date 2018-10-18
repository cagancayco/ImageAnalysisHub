function threshold = previewThreshold(app, EventData)
warning off
if app.PooledButton.Value
    nCells = length(EventData);
    
    for i = 1:nCells
        warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
        currCell = EventData(i);
        
        if length(currCell.events) < 3
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
    
    periods = {EventData.IEI}';
    periods = cellfun(@(X) X, periods, 'UniformOutput', 0);
    periods = cell2mat(periods(~cellfun('isempty', periods)));
    
    if app.TruncatedGaussianSingleButton.Value
        threshold = thTruncGaussianSingle(periods, 0);
    elseif app.logTruncatedGaussianSingleButton.Value
        threshold = thTruncGaussianSingle(log(periods), 1);
    elseif app.GaussianOnlyIntersectionButton.Value
        [~, threshold] = thGaussianOnly(periods, 0);
    elseif app.logGaussianOnlyIntersectionButton.Value
        [~, threshold] = thGaussianOnly(log(periods), 1);
    elseif app.GaussianOnlyMinimumButton.Value
        [threshold, ~] = thGaussianOnly(periods, 0);
    elseif app.logGaussianOnlyMinimumButton.Value
        [threshold, ~] = thGaussianOnly(log(periods), 1);  
    elseif app.GaussianExponentialIntersectionButton.Value
        [~, threshold] = thGaussExp(periods);
    elseif app.logGaussianExpExponentialIntersectionButton.Value
        [~, threshold] = thGaussExpExp(log(periods));
    elseif app.GaussianExponentialMinimumButton.Value
        [threshold, ~] = thGaussExp(periods);
    elseif app.logGaussianExpExponentialMinimumButton.Value
        [threshold, ~] = thGaussExpExp(log(periods));    
    elseif app.BodovaButton.Value
        threshold = thBodova(periods);
    end
else
    threshold = 0;   
end




end