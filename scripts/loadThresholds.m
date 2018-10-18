function loadThresholds(app)
    load(fullfile(app.ThresholdFilePath,app.LoadThresholdResultsEditField.Value));
    
    app.ThresholdResultsTable.Data = cell2table(thresholdsCellArray);
    app.ThresholdResultsTable.ColumnName = varNames;







end