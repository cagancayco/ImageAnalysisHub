function data_struct = loadThresholdData(app)


    inputDir  = app.InputDirectoryEditField.Value;
    slices    = app.AddedSlicesListBox.Items;
    alg_minEASE = app.minEASEButton.Value;
    
    frameRate = app.FrameRateHzEditField.Value;
    startTime = app.StartTimesEditField.Value;
    endTime   = app.EndTimesEditField.Value;
    classesToInc = app.AddedClassesListBox.Items;
    
    data_struct = struct;
    counter = 1;
    
    for i = 1:length(slices)
        
        sliceDir = dir(fullfile(inputDir, slices{i}));
        if alg_minEASE
            try
                for j = 1:length(sliceDir)
                    if ~strcmp(sliceDir(j).name(1),'.')
                        data_struct(counter).sliceName = slices{i};
                        data_struct(counter).cellName  = sliceDir(j).name;
                    
                        fileName = dir(fullfile(inputDir,slices{i},sliceDir(j).name,'*output.mat'));
                        fileName = fullfile(inputDir,slices{i},sliceDir(j).name,fileName.name);
                        data_struct(counter).eventFile = fileName;
                        
                        traceFile = fullfile(inputDir,slices{i},sliceDir(j).name,[sliceDir(j).name,'.mat']);
                        load(traceFile);
                        data = load(fileName);
                        data_struct(counter).algorithm = 'minEASE';
                        if ~isempty(data.eventInfo)
                            samplingInterval = 1/frameRate;
                            data_struct(counter).time = 0:samplingInterval:((data.nSamples-1) * samplingInterval);
                            data_struct(counter).eventClass = data.eventClass;
                            data_struct(counter).events  = data.eventInfo(:,2);
                            data_struct(counter).nSamples   = data.nSamples;
                            data_struct(counter).trace = dataThisCell;
                            data_struct(counter).frameRate = frameRate;
                            data_struct(counter).startTime = startTime;
                            data_struct(counter).endTime = endTime;
                            
                            
                            
                            data_struct(counter).classesToInc = [];
                            for k = 1:length(classesToInc)
                                data_struct(counter).classesToInc = [data_struct(counter).classesToInc,str2double(classesToInc{k})];
                            end
                        end
                        
                        
                        counter = counter + 1;
                    end
                    
                end
            catch
                msgbox('Wrong detection algorithm selected.', 'Incorrect Detection Algorithm', 'error')
                writeEventData = 0;
                break
            end
            
            
        % letsDetect    
        else
            
            try
                fileName = dir(fullfile(inputDir,slices{i},'*ORAMA.mat'));
                fileName = fullfile(inputDir,slices{i},fileName.name);
            
            
                data = load(fileName);
            
            
                for j = 1:length(data.detectoramaOUT.events)
                    data_struct(counter).sliceName = slices{i};
                    data_struct(counter).cellName = [slices{i},'_cell',num2str(j)];
                    data_struct(counter).eventFile = fileName;
            
                    data_struct(counter).algorithm = 'letsDetect';
                    if ~isempty(data.detectoramaOUT.events{j})
                        data_struct(counter).time      = data.detectoramaOUT.timeC;
                        data_struct(counter).events    = unique(sort(data.detectoramaOUT.events{j}(:,1)));
                        data_struct(counter).nSamples  = length(data_struct(counter).events);
                        data_struct(counter).trace     = data.detectoramaOUT.traces(:,j);
                        data_struct(counter).frameRate = frameRate;
                        data_struct(counter).startTime = startTime;
                        data_struct(counter).endTime = endTime;
                           
                    end
                
                    counter = counter + 1;
                end
            
            catch
                %msgbox('Wrong detection algorithm selected.', 'Incorrect Detection Algorithm', 'error')
                break
            end
            
        end

    end

    [~, sortedIdx] = natsort({data_struct.cellName});
    data_struct = data_struct(sortedIdx);
    
end