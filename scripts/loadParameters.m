function loadParameters(app)
    load(fullfile(app.FilePath,app.LoadEventDataMATfileEditField.Value));

    datasets = unique({EventData.sliceName});
    addedSlicesText = '';
    for i = 1:length(datasets)
        addedSlicesText = [addedSlicesText,datasets{i},newline];    
    end

    app.AddedSlicesTextArea.Value = addedSlicesText;
    
    app.DetectionMethodTextArea.Value = EventData(1).algorithm;
    
    if strcmp(EventData(1).algorithm,'minEASE')
        for i = 1:length(EventData(1).classesToInc)
            if i == 1
                app.EventClassesTextArea.Value = num2str(EventData(1).classesToInc(i));
            else
                app.EventClassesTextArea.Value = strcat(app.EventClassesTextArea.Value,', ',num2str(EventData(1).classesToInc(i)));
            end
        end
    end
    
    app.BinTimesTextArea.Value = num2str(EventData(1).binTime);
    
    app.StartTimesTextArea.Value = num2str(EventData(1).startTime);
    
    app.EndTimesTextArea.Value = num2str(EventData(1).endTime);

    app.AnalysisWindowsTable.Data = EventData(1).analysisWindows;
    
    
end