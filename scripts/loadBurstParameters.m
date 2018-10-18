function loadBurstParameters(app)

    load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));

    datasets = unique({EventData.sliceName});
    addedSlicesText = '';
    for i = 1:length(datasets)
        addedSlicesText = [addedSlicesText,datasets{i},newline];    
    end

    app.AddedSlicesTextArea.Value = addedSlicesText;
    
    app.DetectionMethodEditField.Value = EventData(1).algorithm;
    
    if strcmp(EventData(1).algorithm,'minEASE')
        for i = 1:length(EventData(1).classesToInc)
            if i == 1
                app.EventClassesEditField.Value = num2str(EventData(1).classesToInc(i));
            else
                app.EventClassesEditField.Value = strcat(app.EventClassesEditField.Value,', ',num2str(EventData(1).classesToInc(i)));
            end
        end
    end
    
    app.StartTimesEditField.Value = num2str(EventData(1).startTime);
    
    app.EndTimesEditField.Value = num2str(EventData(1).endTime);

    app.AnalysisWindowsTable.Data = EventData(1).analysisWindows;

    app.ThresholdMethodEditField.Value = EventData(1).thresholdMethod;

end