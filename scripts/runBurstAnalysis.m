function runBurstAnalysis(app)
    warning off;
    clc
    if ~ispc
        xlwrite_java_path = [pwd, '/poi_library'];
        javaaddpath(xlwrite_java_path)
    end
    warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
    
    %%% Waitbar stuff %%%
    f = waitbar(0, 'Start burst analysis...', 'Name', 'Burst Analysis Progress');
    
    
    % Create Output filename timestamp appendix
    now = string(datetime('now','Format','d_MMM_y_HH.mm'));
    matFilename = strcat(app.LoadEventDataMATfileEditField.Value(1:end-4),'_',now,'.mat');
    excelFilename = strcat('_BurstAnalysis_',now,'.xlsx');
    
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
    
    newEventData = runSingleThreshold(app, newEventData);

    for i = 1:nCells
        warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
        
        currCell = newEventData(i);
        newEventData(i).burstThreshold = app.MinimumEventsBurstEditField.Value;
        %%% Waitbar stuff %%%
        waitbar(i/(nCells+numel(slices)),f,strrep(sprintf('Analyzing %s', newEventData(i).cellName), '_', '\_'));
        
        
        if length(newEventData(i).events) < 3
            runCell = 0;
        else
            runCell = 1;
        end

        
        switch runCell
        
            case 1
                if app.ManualButton.Value
                    newEventData(i).thresholdMethod = 'Manual';
                    newEventData(i).threshold = app.ManualEditField.Value;
                end
                
                newEventData(i).Bursts = FindBursts(newEventData(i));
                if app.NumberofBurstsCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                        newEventData(i).nBursts = NumBursts(newEventData(i));
                    else
                        newEventData(i).nBursts = [];
                    end
                end
                
                if app.EventsBurstCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                        newEventData(i).nEventsPerBurst = EventsPerBurst(newEventData(i));
                    else
                        newEventData(i).nEventsPerBurst = [];
                    end
                end
                
                if app.BurstDurationCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                    	newEventData(i).burstDuration = BurstDuration(newEventData(i));
                    else
                        newEventData(i).burstDuration = [];
                    end
                end
                
                if app.BurstFrequencyCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                        newEventData(i).burstFreq = BurstFrequency(newEventData(i));
                    else
                        newEventData(i).burstFreq = [];
                    end          
                end
                
                if app.InterburstIntervalCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                        newEventData(i).interburstInt = InterBurstInterval(newEventData(i));
                    else
                        newEventData(i).interburstInt = [];
                    end 
                end
                
                if app.FractionofSingleEventsCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                        newEventData(i).fracSingleEvents = FractionSingleEvents(newEventData(i));
                    else
                        if ~isempty(EventTimes(newEventData(i)))
                            newEventData(i).fracSingleEvents = 1;
                        else
                            newEventData(i).fracSingleEvents = [];
                        end
                    end
                end
                
                if app.TimetoFirstBurstCheckBox.Value
                    if ~isempty(newEventData(i).Bursts)
                        newEventData(i).timeFirstBurst = TimeFirstBurst(newEventData(i));
                    else
                        newEventData(i).timeFirstBurst = [];
                    end
                end
            case 0
        end
        

    end

    EventData = newEventData;
    save(fullfile(app.OutputDirectoryEditField.Value,'MAT-files',matFilename), 'EventData');

    % Create spreadsheet, write summary stats

    % Parameters/Metadata

    parametersHeaders = {'Parameters'; 'Frame Rate (Hz)'; 
                         'Start Time (s)'; 'End Time (s)';
                         'Event Classes'; 'Minimum Events/Burst';
                         'Threshold'; EventData(1).thresholdMethod};

    parametersValues = cell(8,100);
    parametersValues{2,1} = EventData(1).frameRate;
    parametersValues{3,1} = EventData(1).startTime;
    parametersValues{4,1} = EventData(1).endTime;

    if strcmp(EventData(1).algorithm, 'minEASE')
        parametersValues{5,1} = num2str(EventData(1).classesToInc);
    end

    parametersValues{6,1} = EventData(1).burstThreshold;

    for j = 1:numel(slices)
        
        %%% Waitbar stuff %%%
        waitbar((i+j)/(nCells+numel(slices)),f,strrep(sprintf('Writing %s results to Excel', slices{j}), '_', '\_'))
        currSlice = slices{j};
        sliceCells = find(strcmp({EventData.sliceName},currSlice));
        
        sliceFile = strcat(app.OutputDirectoryEditField.Value,'/',currSlice,excelFilename);

        for k = 1:numel(sliceCells)
            parametersValues{7,k} = EventData(sliceCells(k)).cellName(length(currSlice)+2:end);
            parametersValues{8,k} = EventData(sliceCells(k)).threshold;
        end
        
        parameters = horzcat(parametersHeaders, parametersValues);
        
        if ispc
            xlswrite(sliceFile, parameters, 'Parameters');
        else
            xlwrite(sliceFile, parameters, 'Parameters');
        end
        
        
        %% Number of Bursts
        if app.NumberofBurstsCheckBox.Value
            NumBurstsXL(EventData, currSlice, sliceCells, sliceFile)
        end
        %% Events Per Burst
        if app.EventsBurstCheckBox.Value
            EventsBurstXL(EventData, currSlice, sliceCells, sliceFile)
        end
        %% Burst Duration
        if app.BurstDurationCheckBox.Value
            BurstDurationXL(EventData, currSlice, sliceCells, sliceFile)
        end
        %% Burst Frequency
        % Intraburst Frequency
        if app.BurstFrequencyCheckBox.Value
            IntraBurstFreqXL(EventData, currSlice, sliceCells, sliceFile)
        end
        %% Interburst Interval
        if app.InterburstIntervalCheckBox.Value
            InterBurstIntXL(EventData, currSlice, sliceCells, sliceFile)
        end
        %% Fraction of Single Events
        if app.FractionofSingleEventsCheckBox.Value
            FractionSingleEventsXL(EventData, currSlice, sliceCells, sliceFile)
        end
        
        %% Time to First Burst
        if app.TimetoFirstBurstCheckBox.Value
            TimeFirstBurstXL(EventData, currSlice, sliceCells, sliceFile)
        end
    end
    
    delete(f)
end