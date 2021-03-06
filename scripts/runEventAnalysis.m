function runEventAnalysis(app)
    warning off;
    clc
    if ~ispc
        xlwrite_java_path = [pwd, '/poi_library'];
        javaaddpath(xlwrite_java_path)
    end
    warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
    
    %%% Waitbar stuff %%%
    f = waitbar(0, 'Start event analysis...', 'Name', 'Event Analysis Progress');
    
    
    % Create Output filename timestamp appendix
    now = string(datetime('now','Format','d_MMM_y_HH.mm'));
    matFilename = strcat(app.LoadEventDataMATfileEditField.Value(1:end-4),'_',now,'.mat');
    excelFilename = strcat('_EventAnalysis_',now,'.xlsx');


    % Load Event Data MAT-file
    EventData = load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));
    nCells = length(EventData.EventData);
    
    newEventData = EventData.EventData;
    slices = unique({newEventData.sliceName});
    
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
    
    % Load Threshold Results files
    %if ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual') && ~isempty(app.LoadThresholdResultsEditField.Value)
    %    load(fullfile(app.ThresholdFilePath,app.LoadThresholdResultsEditField.Value));
    %    thresholds = cell2table(thresholdsCellArray, 'VariableNames', varNames);
    %elseif ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual') && isempty(app.LoadThresholdResultsEditField.Value)
        newEventData = runSingleThreshold(app, newEventData);
    %end
    
    
    
    
    for i = 1:nCells
        warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
        
        currCell = newEventData(i);
        %%% Waitbar stuff %%%
        waitbar(i/(nCells+numel(slices)),f,strrep(sprintf('Analyzing %s', newEventData(i).cellName), '_', '\_'));
        
        if app.IncludeinAveragesButton.Value
            if length(newEventData(i).events) < 3
                runCell = 0;
            else
                runCell = 1;
            end
        elseif app.RemovefromAveragesButton.Value
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
        end
        
        switch runCell
        
            case 1

                if ~isfield(newEventData(i),'threshold') || isempty(newEventData(i).threshold)
                    % Determine Gap Threshold
                    if app.PooledButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
                        thresholdRow = thresholds(1,:);
                        thresholdGroup = 'Pooled_';
                    elseif app.ConditionLevelButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
                        dataSourceLabels = lower(thresholds.dataSourceLabels);
                        cellRow = find(contains(dataSourceLabels, lower(newEventData(i).sliceName(1:7))));
                        thresholdRow = thresholds(cellRow(1),:);
                        thresholdGroup = 'Condition_';
                    elseif app.ExperimentLevelButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
                        dataSourceLabels = lower(thresholds.dataSourceLabels);
                        cellRow = find(contains(dataSourceLabels, lower(newEventData(i).sliceName)));
                        thresholdRow = thresholds(cellRow(1),:);
                        thresholdGroup = 'Experiment_';
                    elseif app.CellLevelButton.Value && ~strcmp(app.ThresholdingMethodsButtonGroup.SelectedObject.Text,'Manual')
                        dataSourceLabels = thresholds.dataSourceLabels;
                        cellRow = find(strcmpi(dataSourceLabels, newEventData(i).cellName));
                        thresholdRow = thresholds(cellRow(1),:);
                        thresholdGroup = 'Cell_';
                    else
                        thresholdGroup = '';
                    end

                    if app.TruncatedGaussianSingleButton.Value
                        threshold = thresholdRow.IEITruncatedGaussSingleThreshold;
                        thresholdMethod = [thresholdGroup,'TruncatedGaussianSingle'];
                    elseif app.logTruncatedGaussianSingleButton.Value
                        threshold = thresholdRow.logIEITruncatedGaussSingleThreshold;
                        thresholdMethod = [thresholdGroup,'logTrucatedGaussianSingle'];
                    elseif app.GaussianOnlyIntersectionButton.Value
                        threshold = thresholdRow.IEIGaussOnlyThresholdIntersection;
                        thresholdMethod = [thresholdGroup,'GaussianOnlyIntersection'];
                    elseif app.logGaussianOnlyIntersectionButton.Value
                        threshold = thresholdRow.logIEIGaussOnlyThresholdIntersection;
                        thresholdMethod = [thresholdGroup,'logGaussianOnlyIntersection'];
                    elseif app.GaussianOnlyMinimumButton.Value
                        threshold = thresholdRow.IEIGaussOnlyThresholdMinimum;
                        thresholdMethod = [thresholdGroup,'GaussianOnlyMinimum'];
                    elseif app.logGaussianOnlyMinimumButton.Value
                        threshold = thresholdRow.logIEIGaussOnlyThresholdMinimum;  
                        thresholdMethod = [thresholdGroup,'logGaussianOnlyMinimum'];
                    elseif app.GaussianExponentialIntersectionButton.Value
                        threshold = thresholdRow.IEIGaussExpThresholdIntersection;
                        thresholdMethod = [thresholdGroup,'GaussianExponentialIntersection'];
                    elseif app.logGaussianExpExponentialIntersectionButton.Value
                        threshold = thresholdRow.logIEIGaussExpExpThresholdIntersection;
                        thresholdMethod = [thresholdGroup,'logGaussuanExpExponentialIntersection'];
                    elseif app.GaussianExponentialMinimumButton.Value
                        threshold = thresholdRow.IEIGaussExpThresholdIntersection;
                        thresholdMethod = [thresholdGroup,'GaussianExponentialMinimum'];
                    elseif app.logGaussianExpExponentialMinimumButton.Value
                        threshold = thresholdRow.logIEIGaussExpExpThresholdIntersection;    
                        thresholdMethod = [thresholdGroup,'logGaussianExpExponentialMinimum'];
                    elseif app.BodovaButton.Value
                        threshold = thresholdRow.IEIBodovaThreshold;
                        thresholdMethod = [thresholdGroup,'Bodova'];
                    elseif app.ManualButton.Value
                        threshold = app.ManualEditField.Value;
                        thresholdMethod = 'Manual';
                    end

                    newEventData(i).thresholdMethod = thresholdMethod;
                    newEventData(i).threshold = threshold;
                end
                
                % Run event analysis on cell level
                if app.EventsBinCheckBox.Value
                    newEventData(i).EventsPerBin = EventsPerBin(newEventData(i));
                end
        
                if app.FrequencyBinCheckBox.Value
                    newEventData(i).FrequencyPerBin = FrequencyPerBin(newEventData(i));
                end
        
                if app.EventTimingCheckBox.Value || app.CreateclusterFuncOutputCheckBox.Value
                    newEventData(i).EventTiming = EventTimes(newEventData(i));
                end
        
                if app.AllPeriodsCheckBox.Value
                    newEventData(i).AllPeriods = AllPeriods(newEventData(i));
                    
                    if ~isempty(newEventData(i).analysisWindows)
                        newEventData(i).AllPeriods_Windows = AllPeriods_Windows(newEventData(i));
                    end
                end
        
                if app.EventPeriodsCheckBox.Value
                    newEventData(i).EventPeriods = EventPeriods(newEventData(i), newEventData(i).threshold);
                    
                    if ~isempty(newEventData(i).analysisWindows)
                        newEventData(i).EventPeriods_Windows = EventPeriods_Windows(newEventData(i), newEventData(i).threshold);
                    end
                end
        
                if app.QuiescentPeriodsCheckBox.Value
                    newEventData(i).QuiescentPeriods = QuiescentPeriods(newEventData(i), newEventData(i).threshold);
                    
                    if ~isempty(newEventData(i).analysisWindows)
                        newEventData(i).QuiescentPeriods_Windows = QuiescentPeriods_Windows(newEventData(i), newEventData(i).threshold);
                    end
                end
        
                if app.FractionalActiveTimeTimeWindowsCheckBox.Value
                    newEventData(i).FracATWindows = FracATWindows(newEventData(i), newEventData(i).threshold);
                end
        
                if app.FractionalActiveTimeFirstLastEventCheckBox.Value
                    newEventData(i).FracATFirstLast = FracATFirstLast(newEventData(i), newEventData(i).threshold);
                end
                
                %disp([currCell.cellName,': It worked; ',num2str(threshold)])
               

            case 0
                %disp([currCell.cellName,': No events'])
        end
        

    end

    EventData = newEventData;
    save(fullfile(app.OutputDirectoryEditField.Value,'MAT-files',matFilename), 'EventData');


    
    
    % Creating spreadsheet, writing summary stats  
    
    % Parameters/Metadata
    
    parametersHeaders = {'Parameters';'Frame Rate (Hz)'; 'Bin Time (s)';
                             'Start Time (s)'; 'End Time (s)';
                             'Window 1'; 'Window 2'; 'Event Classes'; 'Threshold';
                             EventData(1).thresholdMethod};
                         
    parametersValues = cell(10,100);
    parametersValues{2,1} = EventData(1).frameRate;
    parametersValues{3,1} = EventData(1).binTime;
    parametersValues{4,1} = EventData(1).startTime;
    parametersValues{5,1} = EventData(1).endTime;
    if ~isempty(EventData(1).analysisWindows)
        parametersValues{6,1} = EventData(1).analysisWindows.Var2(1);
        parametersValues{6,2} = EventData(1).analysisWindows.Var3(1);
        parametersValues{7,1} = EventData(1).analysisWindows.Var2(2);
        parametersValues{7,2} = EventData(1).analysisWindows.Var3(2);
    end
    
    if strcmp(EventData(1).algorithm, 'minEASE')
        parametersValues{8,1} = num2str(EventData(1).classesToInc);
    end
    
    
    
    
    
    % Writing parameters to Master Excel File
    %[~,~,master_list_output] = xlsread(app.MasterExcelFileEditField.Value,'Parameters','A:L');
    
    
    for j = 1:numel(slices)
        
        %%% Waitbar stuff %%%
        waitbar((i+j)/(nCells+numel(slices)),f,strrep(sprintf('Writing %s results to Excel', slices{j}), '_', '\_'))
        
        currSlice = slices{j};
        sliceCells = find(strcmp({EventData.sliceName},currSlice));
        
        sliceFile  = strcat(app.OutputDirectoryEditField.Value,'/',currSlice,excelFilename);
        
        for k = 1:numel(sliceCells)
            parametersValues{9,k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            parametersValues{10,k} = EventData(sliceCells(k)).threshold;
        end
        
        parameters = horzcat(parametersHeaders, parametersValues);
        
%         master_list_output_new = cell(1, size(master_list_output,2));
%         master_list_output_new{1,1} = sliceFile;
%         master_list_output_new{1,2} = currSlice;
%         master_list_output_new{1,3} = EventData(1).frameRate;
%         master_list_output_new{1,4} = EventData(1).binTime;
%         master_list_output_new{1,5} = EventData(1).startTime;
%         master_list_output_new{1,6} = EventData(1).endTime;
%         if ~isempty(EventData(1).analysisWindows)
%             master_list_output_new{1,7} = EventData(1).analysisWindows.Var2(1);
%             master_list_output_new{1,8} = EventData(1).analysisWindows.Var3(1);
%             master_list_output_new{1,9} = EventData(1).analysisWindows.Var2(2);
%             master_list_output_new{1,10} = EventData(1).analysisWindows.Var3(2);
%         end
%         if strcmp(EventData(1).algorithm, 'minEASE')
%             master_list_output_new{1,11} = num2str(EventData(1).classesToInc);
%         end
%         master_list_output_new{1,12} = EventData(1).thresholdMethod;
%         
%         master_list_output = [master_list_output; master_list_output_new];
        
        if ispc
            xlswrite(sliceFile, parameters, 'Parameters');
        else
            xlwrite(sliceFile, parameters, 'Parameters');
        end
        
        %% Events Per Bin
        if app.EventsBinCheckBox.Value
            EventsBinXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% Frequency Per Bin
        if app.FrequencyBinCheckBox.Value
            FrequencyBinXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% Event Timing
        if app.EventTimingCheckBox.Value
            EventTimingXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% All Periods
        if app.AllPeriodsCheckBox.Value
            AllPeriodsXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% Event Periods
        if app.EventPeriodsCheckBox.Value
            EventPeriodsXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% Quiescent Periods
        if app.QuiescentPeriodsCheckBox.Value
            QuiescentPeriodsXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% Fractional Active-Time Time Windows
        if app.FractionalActiveTimeTimeWindowsCheckBox.Value
            
            FractionalATFullXL(EventData, currSlice, sliceCells, sliceFile)
            
            if ~isempty(EventData(1).analysisWindows)
                FractionalATWindowsXL(EventData, currSlice, sliceCells, sliceFile)
            end
        end
        
        %% Fractional Active-Time First-Last Event
        if app.FractionalActiveTimeFirstLastEventCheckBox.Value
            FractionalATFirstLastXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% clusterFunc Output
        if app.CreateclusterFuncOutputCheckBox.Value
            CreateClusterFuncOutput(EventData, currSlice, sliceCells, app.OutputDirectoryEditField.Value)
        end
    end
    
%     if ispc
%         xlswrite(app.MasterExcelFileEditField.Value, master_list_output, 'Parameters');
%     else
%         xlwrite(app.MasterExcelFileEditField.Value, master_list_output, 'Parameters');
%     end
    delete(f)
end